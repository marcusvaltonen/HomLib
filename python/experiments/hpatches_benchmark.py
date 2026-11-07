"""
HPatches benchmark of HomLib's point solvers, on precomputed matches (hpatches_features.py).

For each pair, the script runs a robust estimator (LOMSAC by default) with every point solver, and
reports, for the viewpoint and the illumination sequences:
- the share of the true inliers that the estimator finds, the true inliers being the matches within
  the inlier threshold of the ground-truth homography. This share ignores false inliers: the corner
  error accounts for them;
- the AUC of the corner error, which maps the image corners with the estimated model (distortion
  included) and with the ground-truth homography. Wadenbäck et al. (3DV 2026) report the AUC at
  1, 3 and 5 px on HPatches, as "mAA (AUC)", probably in pixels of images resized to 480 px
  (shorter side), as glue-factory's evaluation, while the matches files are in pixels of the
  original images: the numbers are then not directly comparable;
- the median time of a run.

Pairs with no true inlier are left out, since no solver can find any.

``--backend gcransac|magsac|superansac`` replaces LOMSAC with another robust estimator, if HomLib was
built with it (see BACKENDS for the settings that it gets).

``--iterations N`` runs the estimator for exactly N iterations, instead of stopping when it is
confident. Several runs give the accuracy as a function of the time, as in Fig. 6 of Wadenbäck et al.

``--output`` saves the settings and the results of every pair to a JSON file, and ``--compare
a.json b.json`` prints the differences between two such runs (e.g. two versions of HomLib), with their
standard errors over the pairs, and warns about the settings that differ.

Usage
-----
    python hpatches_features.py --data_dir data --output data/hpatches_sift.npz
    python hpatches_benchmark.py data/hpatches_sift.npz [--solvers Wadenback3DV2026TwoSided DLT]
                                 [--backend lomsac] [--seeds 3] [--output results.json]
    python hpatches_benchmark.py --compare before.json after.json
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np
from tqdm import tqdm

import homlib
from hpatches_data import load_matches


def true_inliers(pair, threshold):
    """The matches within the threshold of the ground-truth homography."""
    x, y, H = pair["x"], pair["y"], pair["H"]
    p = H @ np.vstack([x, np.ones(x.shape[1])])
    return np.flatnonzero(np.linalg.norm(p[:2] / p[2] - y, axis=0) < threshold)


def map_points(pose, p):
    """The points p of image 1 mapped to image 2 with an estimated model: undistorted with k1,
    mapped with H, distorted with k2 (division model)."""
    k1, k2 = pose.distortion_parameter, pose.distortion_parameter2
    q = np.array(pose.homography) @ np.vstack([p, 1 + k1 * np.sum(p**2, 0)])
    u = q[:2] / q[2]
    disc = 1 - 4 * k2 * np.sum(u**2, 0)
    with np.errstate(invalid="ignore"):
        return np.where(disc >= 0, u * 2 / (1 + np.sqrt(np.maximum(disc, 0))), np.nan)


def corner_error(pose, pair):
    """The mean distance between the image corners mapped by the estimated model and by the
    ground truth, in pixels."""
    # The points are centered on the image center, ((w - 1) / 2, (h - 1) / 2) with OpenCV's pixel
    # centers at integer coordinates: the outer corners of the image are at (±w / 2, ±h / 2).
    w, h = pair["width"] / 2, pair["height"] / 2
    corners = np.array([[-w, w, w, -w], [-h, -h, h, h]])
    g = pair["H"] @ np.vstack([corners, np.ones(4)])
    err = np.mean(np.linalg.norm(map_points(pose, corners) - g[:2] / g[2], axis=0))
    return err if np.isfinite(err) else np.inf


# The robust estimators: their names in HomLib, and the names of their inlier threshold and of their
# iteration bounds. MAGSAC++ has no threshold, but a largest noise level, sigma_max, which also
# selects its final inliers (the matches within sigma_max). Only LOMSAC takes a random seed: the
# others draw theirs at random, so their runs are not reproducible.
BACKENDS = {
    "lomsac": ("LOMSAC", "inlier_threshold", "min_num_iterations", "max_num_iterations"),
    "gcransac": ("GCRANSAC", "threshold", "min_iterations", "max_iterations"),
    "magsac": ("MAGSAC", "sigma_max", "min_iterations", "max_iterations"),
    "superansac": ("SupeRANSAC", "inlier_threshold", "min_iterations", "max_iterations"),
}


def make_estimator(backend, solver, threshold, iterations, seed):
    """The robust estimator, with the inlier threshold, and exactly this many iterations if
    iterations isn't None"""
    cls, threshold_name, min_name, max_name = BACKENDS[backend]
    options = {threshold_name: threshold}
    if iterations is not None:
        options.update({min_name: iterations, max_name: iterations})
    if backend == "lomsac":
        options["random_seed"] = seed
    return getattr(homlib, cls)(solver, **options)


def inlier_indices(stats):
    """The indices of the inliers that the estimator found"""
    return stats.inlier_indices if hasattr(stats, "inlier_indices") else stats.inliers


def point_solvers():
    return [s for s in homlib.SolverId.__members__.values()
            if homlib.Solver(s).feature_type == homlib.FeatureType.Points]


def main():
    parser = argparse.ArgumentParser(description="HPatches benchmark of HomLib's point solvers")
    parser.add_argument("matches", nargs="?", help="The matches file (hpatches_features.py)")
    parser.add_argument("--threshold", type=float, default=3.0, help="Inlier threshold in pixels")
    parser.add_argument("--backend", choices=sorted(BACKENDS), default="lomsac",
                        help="Robust estimator (default: lomsac). The others need HomLib built with them")
    parser.add_argument("--seeds", type=int, default=3, help="Runs per pair (random seeds, for LOMSAC)")
    parser.add_argument("--iterations", type=int, default=None,
                        help="Run the estimator for exactly this many iterations")
    parser.add_argument("--solvers", nargs="+", choices=[s.name for s in point_solvers()], default=None,
                        help="SolverId names (default: every point solver)")
    parser.add_argument("--output", type=str, default=None, help="Save the settings and the results to this JSON file")
    parser.add_argument("--compare", nargs=2, metavar=("A", "B"), default=None,
                        help="Print the differences B - A between two saved runs, and exit")
    args = parser.parse_args()
    if args.compare:
        compare(*args.compare)
        return

    if not args.matches:
        parser.error("the matches file is required")
    pairs = load_matches(Path(args.matches))

    solvers = ([homlib.SolverId.__members__[s] for s in args.solvers] if args.solvers
               else point_solvers())
    # results[solver][pair] = per seed: (share of the true inliers found, corner error, time)
    results = {s.name: {} for s in solvers}
    for key, pair in tqdm(pairs.items(), desc="Pairs"):
        gt = set(true_inliers(pair, args.threshold).tolist())
        if not gt:
            continue
        x, y = np.ascontiguousarray(pair["x"]), np.ascontiguousarray(pair["y"])
        for s in solvers:
            runs = []
            for seed in range(1, args.seeds + 1):
                estimator = make_estimator(args.backend, homlib.Solver(s), args.threshold, args.iterations, seed)
                t0 = time.perf_counter()
                pose, stats = estimator.estimate(x, y)
                elapsed = time.perf_counter() - t0
                runs.append((len(gt & set(inlier_indices(stats))) / len(gt),
                             corner_error(pose, pair), elapsed))
            results[s.name][key] = runs
    if args.output:
        settings = dict(matches=Path(args.matches).name, threshold=args.threshold, backend=args.backend,
                        seeds=args.seeds, iterations=args.iterations)
        with open(args.output, "w") as f:
            json.dump(dict(settings=settings, results=results), f)
    print_results(results)


def load_results(path):
    """(settings, results) of a file saved with --output. The files saved before the settings were
    the results alone: their settings are None."""
    with open(path) as f:
        data = json.load(f)
    return (data["settings"], data["results"]) if "settings" in data else (None, data)


def warn_settings(settings_a, settings_b):
    """Prints a warning if two runs have different settings, or if one of them has none."""
    if settings_a is None or settings_b is None:
        print("Warning: a file has no settings (saved by an older version of the script): can't check "
              "that the runs have the same settings")
    elif settings_a != settings_b:
        differ = [f"{k} ({settings_a.get(k)} vs {settings_b.get(k)})"
                  for k in sorted(set(settings_a) | set(settings_b)) if settings_a.get(k) != settings_b.get(k)]
        print("Warning: the runs have different settings: " + ", ".join(differ))


def kind(key):
    return "viewpoint" if key.startswith("v_") else "illumination"


def print_results(results):
    print(f"\n{'Solver':<46s} {'Sequences':<13s} {'Inliers found':>13s}  AUC@1  AUC@3  AUC@5 AUC@10  "
          f"{'Time (ms)':>9s}")
    for name, per_pair in sorted(results.items()):
        for k in ("illumination", "viewpoint"):
            runs = [r for key, rs in per_pair.items() if kind(key) == k for r in rs]
            if not runs:
                continue
            recall, err, t = (np.array(v) for v in zip(*runs))
            # AUC of the cumulative error curve up to each threshold
            aucs = [np.mean(np.maximum(0, 1 - err / th)) * 100 for th in (1, 3, 5, 10)]
            print(f"{name:<46s} {k:<13s} {np.mean(recall) * 100:12.2f}%"
                  + "".join(f" {a:6.2f}" for a in aucs) + f"  {np.median(t) * 1e3:9.2f}")


def compare(path_a, path_b):
    """The differences between two runs, b - a, for the pairs and solvers of both: share of the
    true inliers found and AUC@5, averaged over the seeds of each pair, with their standard errors
    over the pairs; and the ratio of the median times, which has no standard error."""
    settings_a, a = load_results(path_a)
    settings_b, b = load_results(path_b)
    warn_settings(settings_a, settings_b)
    print(f"{'Solver':<46s} {'Sequences':<13s} {'Inliers found':>18s} {'AUC@5':>18s} {'Time':>7s}")
    for name in sorted(set(a) & set(b)):
        for k in ("illumination", "viewpoint"):
            keys = [key for key in a[name] if key in b[name] and kind(key) == k]
            if not keys:
                continue
            d_recall = (per_pair(b, name, keys, 0) - per_pair(a, name, keys, 0)) * 100
            d_auc = (per_pair(b, name, keys, 1) - per_pair(a, name, keys, 1)) * 100
            ratio = np.median(per_pair(b, name, keys, 2)) / np.median(per_pair(a, name, keys, 2))
            print(f"{name:<46s} {k:<13s} {d_recall.mean():+8.2f} ±{standard_error(d_recall):5.2f} pts"
                  f" {d_auc.mean():+8.2f} ±{standard_error(d_auc):5.2f} pts {ratio:6.2f}x")


def per_pair(results, name, keys, i):
    """For each pair, the mean over the seeds of: the share of the true inliers found (i = 0), the
    AUC@5 of the corner error (i = 1), or the time (i = 2)"""
    value = [lambda run: run[0], lambda run: max(0.0, 1 - run[1] / 5), lambda run: run[2]][i]
    return np.array([np.mean([value(run) for run in results[name][key]]) for key in keys])


def standard_error(d):
    """The standard error of the mean of d: nan with fewer than 2 values"""
    return d.std(ddof=1) / np.sqrt(len(d)) if len(d) > 1 else np.nan


if __name__ == "__main__":
    main()
