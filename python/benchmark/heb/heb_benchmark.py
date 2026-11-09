"""
HEB benchmark of HomLib's point solvers: the Homography Estimation Benchmark of Barath, Mishkin,
Polic, Förstner and Matas ("A Large Scale Homography Benchmark", CVPR 2023,
https://github.com/danini/homography-benchmark), with the protocol of its paper.

The dataset gives, for homographies between image pairs of 10 scenes of the 1DSfM photo
collections (a pair can have several, of different planes), the SIFT matches with their SNN ratios
(Lowe's ratio: the descriptor distance to the nearest neighbour over the distance to the second
nearest),
which of them are inliers of the ground-truth homography, and the camera intrinsics and relative
pose of each pair (from COLMAP). For each homography and solver, the script keeps
the matches with an SNN ratio below a threshold, runs a robust estimator, and measures, as the
paper (Section 4):
- the rotation error and the absolute translation error (metres, the estimated translation scaled
  to the true baseline) of the pose decomposed from H with the intrinsics (OpenCV's
  decomposeHomographyMat; the solution minimizing the sum of the rotation error and of the
  direction (both in radians) and absolute (metres) translation errors, as the benchmark);
- the transfer error: the mean distance, in image 2, between the ground-truth inliers of image 1
  mapped by the estimated model (distortion included) and their matches (the benchmark's code calls
  it the reprojection error).
A run with fewer than 4 inliers fails (errors of 1e10). An mAA (mean average accuracy) is, for each
of its error thresholds, the share of the homographies with an error below it, averaged over the
thresholds. The metrics are mAAs per scene, averaged over the scenes: rotation over 1 to 10
degrees, translation over 0.1 to 5 m, transfer error over 1 to 20 px (10 thresholds each, the
transfer ones log2-spaced, as the benchmark's code). The main metric, "mAA (10 deg, 5 m)", is the
average of the rotation and translation mAAs. The time is the mean time of a run of the estimator,
per scene then averaged over the scenes, as the paper's Fig. 6.

The protocol has two steps:
- ``tune``: on the training split (2 scenes), with exactly 1000 iterations, a grid search of the
  inlier threshold and of the SNN threshold of each solver (jointly; an SNN threshold of 1 keeps
  every match), for the best main metric (on a tie, the smallest threshold, then the smallest SNN
  threshold); writes them to a settings file. With ``--backend magsac``, the threshold is
  MAGSAC++'s sigma_max, over the same grid. Its scores aren't results: the benchmark's results
  are those of ``test``;
- ``test``: on the test split (the dataset's 8 test scenes), with these settings, each solver run
  for exactly N iterations for each N of ``--iterations`` (by default 10, 100, 1000 and 10000: the
  range of the paper's text; its Fig. 6 goes to 1000), a time-accuracy curve. ``--snn 1`` tests
  without prefiltering, as the paper's other curves. Without ``--settings``, every solver gets a
  4 px threshold and an SNN threshold of 0.8.
The splits are disjoint, as in the paper: dataset_configuration.yaml also lists the 2 training
scenes as test scenes, and they are left out of the test split (the paper's text says nine test
scenes; the dataset's test/ folder has 8).

Differences from the paper and its code:
- Only the two thresholds are tuned (the paper also tunes parameters specific to each estimator),
  over our grid (the paper gives none), on a subset of the homographies (``--max_pairs``).
- The samplers are uniform (LOMSAC's), so the results compare with the paper's curves without
  PROSAC (its scripts default to PROSAC on the SNN-sorted matches).
- A run without an estimation (fewer than 4 matches after prefiltering) has no time, and is left
  out of the mean time; the benchmark's scripts count it as 0. The runs use all the cores
  (``--jobs``), so their times compare within a run.
- The solvers estimate distortion: the points are centred on each image's centre, HomLib's
  distortion centre; the pose is decomposed from H without the distortion; the transfer error maps
  the points with the distortion, and is a failure if a ground-truth inlier has no image (H maps it
  to infinity, or the distortion has no inverse there). HomLib's estimators also keep only the
  solutions with a distortion in their range (k_range).
- The scenes are large (the test split has 8 scenes, 25 GB), so ``--max_pairs`` evaluates a
  random subset of each scene's homographies (despite its name, it counts homographies, not image
  pairs; seeded: the same ones in every run); the paper uses them all (``--max_pairs 0``: days on 4
  cores, and several GB of memory for ``tune``, which keeps the results of every setting).

Requirements: h5py, opencv-python, tqdm, and the dataset from the benchmark's README, in a
directory with ``test/`` and ``training_and_validation/`` (or ``train/``).

Usage
-----
    python heb_benchmark.py tune /path/to/heb --output settings.json [--max_pairs 250]
    python heb_benchmark.py test /path/to/heb [--settings settings.json] [--iterations 10 100 1000 10000]
                                 [--scene NAME] [--max_pairs 1000] [--backend lomsac] [--output heb.json]
    python heb_benchmark.py compare before.json after.json
"""

import argparse
import itertools
import json
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import cv2
import h5py
import numpy as np
from tqdm import tqdm

import homlib

# common.py is in the parent directory (benchmark/)
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import BACKENDS, inlier_indices, make_estimator, map_points, warn_settings  # noqa: E402

# The scenes of each split and the scale of their reconstruction (to metres), as the benchmark's
# dataset_configuration.yaml
SCENES = {
    "train": {"NYC_Library": 6.95844, "Alamo": 3.97745},
    "test": {"Yorkminster": 13.26471, "Tower_of_London": 15.72399, "Madrid_Metropolis": 14.44161,
             "Ellis_Island": 0.94967, "Roman_Forum": 24.38904, "Vienna_Cathedral": 15.07437,
             "Piazza_del_Popolo": 6.25074, "Union_Square": 6.85125},
}
# The thresholds of the mAAs, as the benchmark's code
ROTATION_THRESHOLDS = np.linspace(1, 10, 10)   # degrees
TRANSLATION_THRESHOLDS = np.linspace(0.1, 5, 10)   # metres
TRANSFER_THRESHOLDS = np.logspace(0, np.log2(20), 10, base=2.0)   # px
FAILED = 1e10   # the errors of a run without a model
# The grid of the tuning
TUNE_THRESHOLDS = [1, 2, 3, 4, 6, 8, 12]   # px
TUNE_SNN = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]   # 1.0: no prefiltering (the paper's best: 0.6 for most methods)
DEFAULT_SETTINGS = dict(threshold=4.0, snn=0.8)


def scene_file(path, split, scene):
    for d in (["test"] if split == "test" else ["training_and_validation", "train"]):
        f = Path(path) / d / f"{scene}_homographies.h5"
        if f.exists():
            return f
    raise FileNotFoundError(f"{scene}_homographies.h5 not found in {path}")


def read_pair(f, pair):
    """The matches (Nx10: x1, y1, x2, y2, angles, scales, SNN ratio, ground-truth inlier), the
    relative pose (3x4) and the intrinsics and sizes (width, height) of the two images"""
    a, b = "_".join(pair.split("_")[0:3]), "_".join(pair.split("_")[3:6])
    return (f[f"corr_{pair}"][()], f[f"pose_{pair}"][()], f[f"K_{a}"][()], f[f"K_{b}"][()],
            np.ravel(f[f"size_{a}"][()]), np.ravel(f[f"size_{b}"][()]))


def pose_errors(H, pose, K1, K2, scale):
    """The rotation error (degrees) and the absolute translation error (metres) of the pose
    decomposed from H (pixels), as the benchmark: the decomposition minimizing the sum of the
    rotation and translation-direction errors (radians) and of the absolute translation error"""
    R_gt, t_gt = pose[:, :3], scale * pose[:, 3]
    _, Rs, ts, _ = cv2.decomposeHomographyMat(np.linalg.inv(K2) @ H @ K1, np.eye(3))
    best, errs = np.inf, (FAILED, FAILED)
    for R, t in zip(Rs, ts):
        t = t.ravel()
        if not (np.all(np.isfinite(R)) and np.all(np.isfinite(t))):
            continue
        err_R = np.arccos(np.clip((np.trace(R_gt @ R.T) - 1) / 2, -1, 1))
        n, n_gt = np.linalg.norm(t), np.linalg.norm(t_gt)
        err_t = np.arccos(min(1.0, abs(t @ t_gt) / (n * n_gt + 1e-15)))   # the direction, up to sign
        err_abs_t = np.linalg.norm(n_gt * t / (n + 1e-15) - t_gt)
        if err_R + err_t + err_abs_t < best:
            best, errs = err_R + err_t + err_abs_t, (np.degrees(err_R), err_abs_t)
    return errs


def run_pair(job):
    """For each configuration (solver, threshold, SNN threshold, iterations): per seed, (rotation
    error, translation error, transfer error, number of inliers, time)"""
    path, split, scene, pair, configs, backend, seeds = job
    with h5py.File(scene_file(path, split, scene), "r") as f:
        corr, pose, K1, K2, size1, size2 = read_pair(f, pair)
    c1, c2 = size1[:2] / 2, size2[:2] / 2   # the sizes are (width, height)
    gt = corr[:, 9].astype(bool)
    T1 = np.array([[1, 0, -c1[0]], [0, 1, -c1[1]], [0, 0, 1.0]])
    T2i = np.array([[1, 0, c2[0]], [0, 1, c2[1]], [0, 0, 1.0]])
    out = {}
    for key in configs:
        name, threshold, snn, iterations = key.split("|")
        keep = corr[:, 8] < float(snn) if float(snn) < 1 else np.ones(len(corr), bool)
        x = np.ascontiguousarray((corr[keep, 0:2] - c1).T)
        y = np.ascontiguousarray((corr[keep, 2:4] - c2).T)
        solver = homlib.Solver(homlib.solver_ids(names=[name])[0])
        runs = []
        for seed in range(1, seeds + 1):
            if x.shape[1] < 4:   # no estimation, so no time (NaN, left out of the mean time)
                runs.append((FAILED, FAILED, FAILED, 0, None))
                continue
            estimator = make_estimator(backend, solver, float(threshold), int(iterations), seed)
            t0 = time.perf_counter()
            model, stats = estimator.estimate(x, y)
            elapsed = time.perf_counter() - t0
            n_in = len(inlier_indices(stats))
            err_R = err_t = transfer = FAILED
            if n_in >= 4:
                err_R, err_t = pose_errors(T2i @ np.array(model.homography) @ T1, pose, K1, K2, SCENES[split][scene])
                if gt.sum() > 1:
                    mapped, ok = map_points(model, (corr[gt, 0:2] - c1).T)
                    if ok.all():
                        transfer = np.mean(np.linalg.norm(mapped - (corr[gt, 2:4] - c2).T, axis=0))
            runs.append((float(err_R), float(err_t), float(transfer), n_in, elapsed))
        out[key] = runs
    return scene, pair, out


def accuracy(errors, thresholds):
    """The share of the thresholds that each error is below, per run"""
    e = np.asarray(errors, dtype=float)[:, None]
    return (e <= thresholds[None, :]).mean(1)


def scores(runs):
    """Per run: the accuracies (main, rotation, translation, transfer), and the times"""
    err_R, err_t, transfer, _, t = (np.array(v, dtype=float) for v in zip(*runs))
    a_R, a_t = accuracy(err_R, ROTATION_THRESHOLDS), accuracy(err_t, TRANSLATION_THRESHOLDS)
    return np.stack([(a_R + a_t) / 2, a_R, a_t, accuracy(transfer, TRANSFER_THRESHOLDS)]), t


def summarize(per_scene):
    """The mAAs (main, rotation, translation, transfer; %) averaged over the scenes, the
    per-scene main mAA, and the mean time of a run, per scene then averaged over the scenes (as
    the paper's Fig. 6)"""
    means, mains, times = [], {}, []
    for scene, per_pair in per_scene.items():
        acc, t = scores([r for rs in per_pair.values() for r in rs])
        means.append(acc.mean(1))
        mains[scene] = 100 * acc[0].mean()
        times.append(np.nanmean(t))
    return 100 * np.mean(means, 0), mains, np.nanmean(times)


def make_jobs(args, split, configs):
    jobs = []
    for scene in ([args.scene] if args.scene else SCENES[split]):
        with h5py.File(scene_file(args.path, split, scene), "r") as f:
            pairs = sorted(k[len("corr_"):] for k in f.keys() if k.startswith("corr_"))
        if 0 < args.max_pairs < len(pairs):
            pairs = sorted(np.random.default_rng(0).choice(pairs, args.max_pairs, replace=False))
        jobs += [(args.path, split, scene, p, configs, args.backend, args.seeds) for p in pairs]
    return jobs


def run(jobs, configs, processes):
    """results[configuration][scene][pair] = per seed: (rotation error, translation error,
    transfer error, inliers, time)"""
    results = {key: {} for key in configs}
    with Pool(processes) as pool:
        for scene, pair, out in tqdm(pool.imap_unordered(run_pair, jobs, chunksize=4), total=len(jobs), desc="Pairs"):
            for key, runs in out.items():
                results[key].setdefault(scene, {})[pair] = runs
    return results


def tune(args):
    configs = [f"{s}|{th}|{snn}|1000" for s, th, snn in itertools.product(args.solvers, TUNE_THRESHOLDS, TUNE_SNN)]
    results = run(make_jobs(args, "train", configs), configs, args.jobs)
    settings = {}
    for name in args.solvers:
        # The best main mAA; on a tie, the first in the order of the grid (the smallest threshold,
        # then the smallest SNN threshold)
        grid = [(summarize(results[f"{name}|{th}|{snn}|1000"])[0][0], th, snn)
                for th, snn in itertools.product(TUNE_THRESHOLDS, TUNE_SNN)]
        best = max(grid, key=lambda g: g[0])
        ties = sum(g[0] == best[0] for g in grid) - 1
        settings[name] = dict(threshold=float(best[1]), snn=float(best[2]), train_mAA=float(best[0]))
        print(f"{name:<36s} threshold {best[1]:>4} px, SNN {best[2]:.2f}: mAA (10 deg, 5 m) {best[0]:.2f}"
              + (f" ({ties} other settings as good)" if ties else ""))
    with open(args.output, "w") as f:
        json.dump(dict(settings=vars_settings(args, "train"), solvers=settings), f, indent=1)


def vars_settings(args, split):
    return dict(split=split, scene=args.scene, max_pairs=args.max_pairs, backend=args.backend, seeds=args.seeds)


def test(args):
    tuned = {}
    if args.settings:
        with open(args.settings) as f:
            tuned = json.load(f)["solvers"]
    missing = [name for name in args.solvers if args.settings and name not in tuned]
    if missing:
        sys.exit(f"{args.settings} has no settings for {', '.join(missing)}: run tune with these solvers")
    configs = []
    for name in args.solvers:
        s = tuned.get(name, DEFAULT_SETTINGS)
        snn = args.snn if args.snn is not None else s['snn']
        configs += [f"{name}|{s['threshold']}|{snn}|{n}" for n in args.iterations]
    results = run(make_jobs(args, "test", configs), configs, args.jobs)
    if args.output:
        settings = dict(vars_settings(args, "test"), settings_file=args.settings)
        with open(args.output, "w") as f:
            json.dump(dict(settings=settings, results=results), f)
    print_results(results, args.per_scene)


def print_results(results, per_scene=False):
    for n in sorted({int(k.split("|")[3]) for k in results}):
        print(f"\n{n} iterations: mAA (%) averaged over the scenes, mean time of a run")
        print(f"{'Solver':<36s} {'Thr.':>5s} {'SNN':>5s} {'(10°, 5 m)':>10s} {'Rot.':>6s} {'Transl.':>7s} "
              f"{'Transfer':>8s} {'Time (ms)':>9s}")
        for key in sorted(k for k in results if int(k.split("|")[3]) == n):
            name, th, snn, _ = key.split("|")
            m, mains, t = summarize(results[key])
            print(f"{name:<36s} {float(th):5.1f} {float(snn):5.2f} {m[0]:10.2f} {m[1]:6.2f} {m[2]:7.2f} {m[3]:7.2f} "
                  f"{t * 1e3:9.2f}")
            if per_scene:
                print("    " + ", ".join(f"{s} {v:.1f}" for s, v in sorted(mains.items())))


def compare(path_a, path_b):
    """For the configurations of both runs (same solver, settings and iterations): the differences
    B - A of the main and transfer mAAs, per scene then averaged over the scenes as the mAAs,
    with their standard errors (from the per-homography differences), and the ratio of the mean
    times"""
    with open(path_a) as f:
        a = json.load(f)
    with open(path_b) as f:
        b = json.load(f)
    warn_settings(a.get("settings"), b.get("settings"))
    ra, rb = a["results"], b["results"]
    other = sorted({k.split("|")[0] for k in ra} & {k.split("|")[0] for k in rb}
                   - {k.split("|")[0] for k in set(ra) & set(rb)})
    if other:
        print("Warning: no configuration (threshold, SNN, iterations) in common in the two runs, not compared: "
              + ", ".join(other))
    print(f"{'Solver':<36s} {'Iter.':>6s} {'mAA (10°, 5 m)':>18s} {'Transfer':>18s} {'Time':>6s}")
    for key in sorted(set(ra) & set(rb), key=lambda k: (k.split("|")[0], int(k.split("|")[3]))):
        means, variances, ta, tb = [], [], [], []
        scenes = sorted(set(ra[key]) & set(rb[key]))
        for scene in scenes:
            d, sa, sb = [], [], []
            for pair in sorted(set(ra[key][scene]) & set(rb[key][scene])):
                (acc_a, t_a), (acc_b, t_b) = scores(ra[key][scene][pair]), scores(rb[key][scene][pair])
                d.append(acc_b.mean(1) - acc_a.mean(1))
                sa += list(t_a)
                sb += list(t_b)
            if not d:   # no homography of this scene in both runs
                continue
            ta.append(np.nanmean(sa))   # the mean times per scene, averaged over the scenes as in summarize
            tb.append(np.nanmean(sb))
            d = 100 * np.array(d)
            means.append(d.mean(0))
            variances.append(d.var(0, ddof=1) / len(d) if len(d) > 1 else np.zeros(4))
        mean, se = np.mean(means, 0), np.sqrt(np.sum(variances, 0)) / len(means)
        name, n = key.split("|")[0], key.split("|")[3]
        print(f"{name:<36s} {n:>6s} {mean[0]:+8.2f} ±{se[0]:5.2f} pts {mean[3]:+8.2f} ±{se[3]:5.2f} pts "
              f"{np.nanmean(tb) / np.nanmean(ta):5.2f}x")


def main():
    parser = argparse.ArgumentParser(description="HEB benchmark of HomLib's point solvers")
    sub = parser.add_subparsers(dest="command", required=True)
    point_solvers = [s.name for s in homlib.solver_ids(homlib.FeatureType.Points)]
    for name in ("tune", "test"):
        p = sub.add_parser(name)
        p.add_argument("path", help="The dataset directory (with test/ and training_and_validation/)")
        p.add_argument("--scene", default=None, help="One scene of the split (default: all)")
        p.add_argument("--backend", choices=sorted(BACKENDS), default="lomsac",
                       help="Robust estimator (default: lomsac). The others need HomLib built with them")
        p.add_argument("--seeds", type=int, default=1, help="Runs per pair (random seeds, for LOMSAC)")
        p.add_argument("--solvers", nargs="+", choices=point_solvers, default=point_solvers)
        p.add_argument("--jobs", type=int, default=os.cpu_count(), help="Processes (default: one per core)")
    tune_p, test_p = sub.choices["tune"], sub.choices["test"]
    tune_p.add_argument("--max_pairs", type=int, default=250,
                        help="Homographies per scene, chosen at random (seeded), 0 for all (default: 250)")
    tune_p.add_argument("--output", required=True, help="The settings file to write")
    test_p.add_argument("--max_pairs", type=int, default=1000,
                        help="Homographies per scene, chosen at random (seeded), 0 for all (default: 1000)")
    test_p.add_argument("--settings", default=None, help="The settings file written by tune")
    test_p.add_argument("--snn", type=float, default=None,
                        help="SNN threshold for every solver, instead of the tuned one (1: no prefiltering, as the "
                             "paper's curves without it)")
    test_p.add_argument("--iterations", type=int, nargs="+", default=[10, 100, 1000, 10000],
                        help="Run each solver for exactly this many iterations, for each value (default: 10 100 1000 "
                             "10000)")
    test_p.add_argument("--per_scene", action="store_true", help="Also print the main mAA of each scene")
    test_p.add_argument("--output", default=None, help="Save the settings and the results to this JSON file")
    cmp_p = sub.add_parser("compare")
    cmp_p.add_argument("files", nargs=2, metavar=("A", "B"), help="Print the differences B - A between two test runs")
    args = parser.parse_args()
    if args.command == "compare":
        compare(*args.files)
        return
    split = "train" if args.command == "tune" else "test"
    if args.scene and args.scene not in SCENES[split]:
        parser.error(f"{args.scene} isn't in the {split} split: {', '.join(SCENES[split])}")
    tune(args) if args.command == "tune" else test(args)


if __name__ == "__main__":
    main()
