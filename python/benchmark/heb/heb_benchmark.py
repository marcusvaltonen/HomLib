"""
HEB benchmark of HomLib's point solvers: the Homography Estimation Benchmark of Barath, Mishkin,
Polic, Förstner and Matas ("A Large Scale Homography Benchmark", CVPR 2023,
https://github.com/danini/homography-benchmark).

The dataset gives, for image pairs of 10 scenes of the 1DSfM photo collections, the SIFT matches
with their SNN ratios, which of them are inliers of the ground-truth homography, and the camera
intrinsics and relative pose of each pair (from COLMAP). Like the benchmark's own scripts, this
one keeps the matches with an SNN ratio below ``--snn`` (0.8), runs a robust estimator on them,
and reports, for each solver:
- mAA of the reprojection error: the mean distance between the ground-truth inliers of image 2 and
  their mapping by the estimated model (distortion included), averaged over the thresholds
  2^linspace(0, log2(20), 10) px, as the benchmark;
- mAA of the pose error: the larger of the rotation and translation (direction) errors of the
  pose decomposed from H (OpenCV's decomposeHomographyMat, with the intrinsics; the best of its
  solutions), over the thresholds 1, 2, ..., 10 degrees, as the benchmark;
- the median time of a run of the estimator.
A run with fewer than 4 inliers gets errors of 1e10, as in the benchmark. The points are centred on each image's
centre, HomLib's distortion centre; the distortion solvers' H is taken without the distortion for
the pose.

The scenes are large (the test split has 8 scenes, 25 GB), so ``--max_pairs`` evaluates a random
subset of each scene's pairs (seeded: the same pairs in every run).

Requirements: h5py, opencv-python, tqdm, and the dataset from the benchmark's README, in a
directory with ``test/`` and ``training_and_validation/`` (or ``train/``).

Usage
-----
    python heb_benchmark.py /path/to/heb [--split test] [--scene NAME] [--max_pairs 1000]
                            [--backend lomsac] [--threshold 4] [--solvers ...] [--output heb.json]
    python heb_benchmark.py --compare before.json after.json
"""

import argparse
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
REPROJECTION_THRESHOLDS = np.logspace(0, np.log2(20), 10, base=2.0)   # px
POSE_THRESHOLDS = np.linspace(1, 10, 10)   # degrees
FAILED = 1e10   # the errors of a run without a model


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


def pose_error(H, pose, K1, K2, scale):
    """The larger of the rotation and translation errors (degrees) of the best pose decomposed from
    H (pixels), as the benchmark: the decomposition minimizing the sum of the rotation and
    translation errors (radians) and of the absolute translation error (metres)."""
    R_gt, t_gt = pose[:, :3], scale * pose[:, 3]
    _, Rs, ts, _ = cv2.decomposeHomographyMat(np.linalg.inv(K2) @ H @ K1, np.eye(3))
    best, err = np.inf, FAILED
    for R, t in zip(Rs, ts):
        t = t.ravel()
        if not (np.all(np.isfinite(R)) and np.all(np.isfinite(t))):
            continue
        err_R = np.arccos(np.clip((np.trace(R_gt @ R.T) - 1) / 2, -1, 1))
        n, n_gt = np.linalg.norm(t), np.linalg.norm(t_gt)
        err_t = np.arccos(min(1.0, abs(t @ t_gt) / (n * n_gt + 1e-15)))   # the direction, up to sign
        err_abs_t = np.linalg.norm(n_gt * t / (n + 1e-15) - t_gt)
        if err_R + err_t + err_abs_t < best:
            best, err = err_R + err_t + err_abs_t, np.degrees(max(err_R, err_t))
    return err


def run_pair(job):
    """Per solver: per seed, (reprojection error, pose error, number of inliers, time)"""
    path, split, scene, pair, args = job
    with h5py.File(scene_file(path, split, scene), "r") as f:
        corr, pose, K1, K2, size1, size2 = read_pair(f, pair)
    c1, c2 = size1[:2] / 2, size2[:2] / 2   # the sizes are (width, height)
    gt = corr[:, 9].astype(bool)
    keep = corr[:, 8] < args.snn
    x = np.ascontiguousarray((corr[keep, 0:2] - c1).T)
    y = np.ascontiguousarray((corr[keep, 2:4] - c2).T)
    T1 = np.array([[1, 0, -c1[0]], [0, 1, -c1[1]], [0, 0, 1.0]])
    T2i = np.array([[1, 0, c2[0]], [0, 1, c2[1]], [0, 0, 1.0]])
    out = {}
    for name in args.solvers:
        runs = []
        for seed in range(1, args.seeds + 1):
            if x.shape[1] < 4:
                runs.append((FAILED, FAILED, 0, 0.0))
                continue
            estimator = make_estimator(args.backend, homlib.Solver(homlib.solver_ids(names=[name])[0]),
                                       args.threshold, args.iterations, seed)
            t0 = time.perf_counter()
            model, stats = estimator.estimate(x, y)
            elapsed = time.perf_counter() - t0
            n_in = len(inlier_indices(stats))
            reproj = perr = FAILED
            if n_in >= 4 and gt.sum() > 1:
                mapped, ok = map_points(model, (corr[gt, 0:2] - c1).T)
                if ok.all():
                    reproj = np.mean(np.linalg.norm(mapped - (corr[gt, 2:4] - c2).T, axis=0))
                perr = pose_error(T2i @ np.array(model.homography) @ T1, pose, K1, K2, SCENES[split][scene])
            runs.append((float(reproj), float(perr), n_in, elapsed))
        out[name] = runs
    return scene, pair, out


def mAA(errors, thresholds):
    """The mean over the pairs of the share of the thresholds that each error is below, in %"""
    e = np.asarray(errors, dtype=float)
    return 100 * np.mean([(e <= th).mean() for th in thresholds]) if e.size else np.nan


def summarize(results):
    """Per solver: (mAA of the reprojection error, mAA of the pose error, mean inliers, median time)"""
    rows = {}
    for name, per_pair in results.items():
        runs = [r for rs in per_pair.values() for r in rs]
        reproj, perr, n_in, t = (np.array(v, dtype=float) for v in zip(*runs))
        rows[name] = (mAA(reproj, REPROJECTION_THRESHOLDS), mAA(perr, POSE_THRESHOLDS), n_in.mean(), np.median(t))
    return rows


def print_results(results):
    n = len(next(iter(results.values()), {}))
    print(f"\n{n} pairs")
    print(f"{'Solver':<36s} {'mAA reproj.':>11s} {'mAA pose':>9s} {'Inliers':>8s} {'Time (ms)':>9s}")
    for name, (ar, ap, n_in, t) in sorted(summarize(results).items()):
        print(f"{name:<36s} {ar:11.2f} {ap:9.2f} {n_in:8.1f} {t * 1e3:9.2f}")


def compare(path_a, path_b):
    """The differences B - A of the two mAAs, with their standard errors over the pairs, and the
    ratio of the median times"""
    with open(path_a) as f:
        a = json.load(f)
    with open(path_b) as f:
        b = json.load(f)
    warn_settings(a.get("settings"), b.get("settings"))
    ra, rb = a["results"], b["results"]
    print(f"{'Solver':<36s} {'mAA reproj.':>18s} {'mAA pose':>18s} {'Time':>6s}")
    for name in sorted(set(ra) & set(rb)):
        keys = sorted(set(ra[name]) & set(rb[name]))
        cols = []
        for i, ths in ((0, REPROJECTION_THRESHOLDS), (1, POSE_THRESHOLDS)):
            d = np.array([np.mean([mAA([r[i] for r in rb[name][k]], ths)]) -
                          np.mean([mAA([r[i] for r in ra[name][k]], ths)]) for k in keys])
            cols.append(f"{d.mean():+8.2f} ±{d.std(ddof=1) / np.sqrt(len(d)) if len(d) > 1 else 0:5.2f} pts")
        ta = np.median([r[3] for k in keys for r in ra[name][k]])
        tb = np.median([r[3] for k in keys for r in rb[name][k]])
        print(f"{name:<36s} {cols[0]:>18s} {cols[1]:>18s} {tb / ta:5.2f}x")


def main():
    parser = argparse.ArgumentParser(description="HEB benchmark of HomLib's point solvers")
    parser.add_argument("path", nargs="?", help="The dataset directory (with test/ and training_and_validation/)")
    parser.add_argument("--split", choices=["test", "train"], default="test")
    parser.add_argument("--scene", choices=sorted(set(SCENES["test"]) | set(SCENES["train"])), default=None,
                        help="One scene (default: every scene of the split)")
    parser.add_argument("--max_pairs", type=int, default=1000,
                        help="Pairs per scene, chosen at random (seeded), 0 for all (default: 1000)")
    parser.add_argument("--snn", type=float, default=0.8, help="SNN ratio threshold of the matches (default: 0.8)")
    parser.add_argument("--threshold", type=float, default=4.0, help="Inlier threshold in pixels (default: 4)")
    parser.add_argument("--backend", choices=sorted(BACKENDS), default="lomsac",
                        help="Robust estimator (default: lomsac). The others need HomLib built with them")
    parser.add_argument("--seeds", type=int, default=1, help="Runs per pair (random seeds, for LOMSAC)")
    parser.add_argument("--iterations", type=int, default=None,
                        help="Run the estimator for exactly this many iterations")
    point_solvers = [s.name for s in homlib.solver_ids(homlib.FeatureType.Points)]
    parser.add_argument("--solvers", nargs="+", choices=point_solvers, default=point_solvers)
    parser.add_argument("--jobs", type=int, default=os.cpu_count(), help="Processes (default: one per core)")
    parser.add_argument("--output", type=str, default=None, help="Save the settings and the results to this JSON file")
    parser.add_argument("--compare", nargs=2, metavar=("A", "B"), default=None,
                        help="Print the differences B - A between two saved runs, and exit")
    args = parser.parse_args()
    if args.compare:
        compare(*args.compare)
        return
    if not args.path:
        parser.error("the dataset directory is required")
    if args.scene and args.scene not in SCENES[args.split]:
        parser.error(f"{args.scene} isn't in the {args.split} split")

    jobs = []
    for scene in ([args.scene] if args.scene else SCENES[args.split]):
        with h5py.File(scene_file(args.path, args.split, scene), "r") as f:
            pairs = sorted(k[len("corr_"):] for k in f.keys() if k.startswith("corr_"))
        if 0 < args.max_pairs < len(pairs):
            rng = np.random.default_rng(0)
            pairs = sorted(rng.choice(pairs, args.max_pairs, replace=False))
        jobs += [(args.path, args.split, scene, p, args) for p in pairs]

    # results[solver][scene/pair] = per seed: (reprojection error, pose error, inliers, time)
    results = {name: {} for name in args.solvers}
    with Pool(args.jobs) as pool:
        for scene, pair, out in tqdm(pool.imap_unordered(run_pair, jobs, chunksize=4), total=len(jobs), desc="Pairs"):
            for name, runs in out.items():
                results[name][f"{scene}/{pair}"] = runs
    if args.output:
        settings = dict(split=args.split, scene=args.scene, max_pairs=args.max_pairs, snn=args.snn,
                        threshold=args.threshold, backend=args.backend, seeds=args.seeds, iterations=args.iterations)
        with open(args.output, "w") as f:
            json.dump(dict(settings=settings, results=results), f)
    print_results(results)


if __name__ == "__main__":
    main()
