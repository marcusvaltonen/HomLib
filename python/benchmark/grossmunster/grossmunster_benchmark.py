"""
Grossmünster Church benchmark of HomLib's two-sided solvers, on the matches and the ground truth of
grossmunster_ground_truth.py: Table 3 of Wadenbäck et al. (3DV 2026).

The images have fisheye distortion, which the one-parameter division model only approximates, in
both images of each pair: the paper evaluates the two-sided solvers only, those with independent
distortions (k1, k2) and those with an equal distortion (k1 = k2). Each solver runs in LOMSAC (by
default), as in the paper: at least 500 iterations (``--min_iterations``), then a final non-linear
least-squares refinement on the inliers. For each solver, the script reports the means over the
pairs and the seeds of:
- the paper's homography error, eps = ||H - H_est||_F, where the ground truth H and the estimate
  H_est are each scaled to unit Frobenius norm and positive determinant, in calibrated coordinates:
  K^-1 H K, where H maps the undistorted points, centered on the principal point, in pixels
  (grossmunster_ground_truth.py), and K is the focal lengths. The paper doesn't say in which
  coordinates it measures eps, which depends on them; its synthetic experiments, as
  synthetic_benchmark.py, measure it for a focal length of 1. The script also gives eps in the
  pixels of the matches ("Error (px)"), where the translation entries of H dominate its norm, so
  that eps is smaller and mostly measures the translation;
- the number of inliers found;
- the time of a run.
``--per_pair`` also prints the errors of each pair (median over the seeds).

``--pairs`` selects the pairs among those of the matches file (default: PAIRS in
grossmunster_data.py, the paper's five known pairs and a guess of its sixth, (2519, 2520); the
other candidate, 2515-2517, is in the file too).

``--backend gcransac|magsac|superansac`` replaces LOMSAC with another robust estimator, if HomLib was
built with it (see BACKENDS in common.py); they run their own final refinement.

``--output`` saves the settings and the results of every pair to a JSON file, and ``--compare
a.json b.json`` prints the differences between two such runs, with their standard errors over the
pairs, and warns about the settings that differ.

Usage
-----
    python grossmunster_ground_truth.py --data_dir data/fisheye_grossmunster --sfm_dir data/grossmunster_sfm \
        --output data/grossmunster_sift.npz
    python grossmunster_benchmark.py data/grossmunster_sift.npz [--threshold 3] [--seeds 10] [--output results.json]
    python grossmunster_benchmark.py --compare before.json after.json
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

import homlib
from grossmunster_data import PAIRS, load_matches, pair_key

# common.py is in the parent directory (benchmark/)
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import BACKENDS, inlier_indices, make_estimator, warn_settings  # noqa: E402

# Table 3 of the paper: the solvers, and their mean error and mean time (ms)
TABLE3 = [("KukelovaCVPR2015TwoSided", "Kukelova et al. (5pt)", 0.039, 2657),
          ("KukelovaCVPR2015TwoSided6Pt", "Kukelova et al. (6pt)", 0.040, 330),
          ("Wadenback3DV2026TwoSided", "Wadenbäck et al. (5pt)", 0.038, 150),
          ("FitzgibbonCVPR2001TwoSidedEqual", "Fitzgibbon (5pt), k1 = k2", 0.040, 233),
          ("KukelovaCVPR2015TwoSidedEqual", "Kukelova et al. (5pt), k1 = k2", 0.040, 2053),
          ("KukelovaCVPR2015TwoSidedEqual6Pt", "Kukelova et al. (6pt), k1 = k2", 0.039, 321),
          ("Wadenback3DV2026TwoSidedEqual", "Wadenbäck et al. (4.5pt), k1 = k2", 0.039, 84)]


def normalized(H):
    """H scaled to unit Frobenius norm and positive determinant"""
    H = H / np.linalg.norm(H)
    return -H if np.linalg.det(H) < 0 else H


def homography_errors(H, H_est, K):
    """The paper's error ||H - H_est||_F of the normalized homographies, in calibrated coordinates
    (K^-1 H K), and in the coordinates of the matches (pixels)"""
    Ki = np.linalg.inv(K)
    return (np.linalg.norm(normalized(Ki @ H @ K) - normalized(Ki @ H_est @ K)),
            np.linalg.norm(normalized(H) - normalized(H_est)))


def main():
    parser = argparse.ArgumentParser(description="Grossmünster Church benchmark of HomLib's two-sided solvers")
    parser.add_argument("matches", nargs="?", help="The matches file (grossmunster_ground_truth.py)")
    parser.add_argument("--threshold", type=float, default=3.0, help="Inlier threshold in pixels (default: 3)")
    parser.add_argument("--backend", choices=sorted(BACKENDS), default="lomsac",
                        help="Robust estimator (default: lomsac). The others need HomLib built with them")
    parser.add_argument("--min_iterations", type=int, default=500, help="Minimum number of iterations (default: 500)")
    parser.add_argument("--seeds", type=int, default=10, help="Runs per pair (random seeds, for LOMSAC)")
    parser.add_argument("--pairs", nargs="+", default=None,
                        help="The pairs, as a-b (default: PAIRS of grossmunster_data.py)")
    parser.add_argument("--solvers", nargs="+", choices=[s for s, *_ in TABLE3], default=None,
                        help="SolverId names (default: the seven of Table 3)")
    parser.add_argument("--per_pair", action="store_true", help="Also print the errors of each pair")
    parser.add_argument("--output", type=str, default=None, help="Save the settings and the results to this JSON file")
    parser.add_argument("--compare", nargs=2, metavar=("A", "B"), default=None,
                        help="Print the differences B - A between two saved runs, and exit")
    args = parser.parse_args()
    if args.compare:
        compare(*args.compare)
        return
    if not args.matches:
        parser.error("the matches file is required")

    pairs, scene = load_matches(Path(args.matches))
    keys = args.pairs or [pair_key(a, b) for a, b in PAIRS]
    missing = [k for k in keys if k not in pairs]
    if missing:
        parser.error(f"pairs not in {args.matches}: {', '.join(missing)} (it has {', '.join(sorted(pairs))})")
    names = args.solvers or [s for s, *_ in TABLE3]
    # The backends other than LOMSAC refine their final model by default
    extra = dict(final_least_squares=True) if args.backend == "lomsac" else {}
    # results[solver][pair] = per seed: (error in calibrated coordinates, error in pixels, inliers, time)
    results = {name: {} for name in names}
    for name in names:
        for key in keys:
            pair = pairs[key]
            x, y = np.ascontiguousarray(pair["x"]), np.ascontiguousarray(pair["y"])
            runs = []
            for seed in range(1, args.seeds + 1):
                estimator = make_estimator(args.backend, homlib.Solver(homlib.SolverId.__members__[name]),
                                           args.threshold, None, seed, min_iterations=args.min_iterations, **extra)
                t0 = time.perf_counter()
                pose, stats = estimator.estimate(x, y)
                elapsed = time.perf_counter() - t0
                runs.append((*homography_errors(pair["H"], np.array(pose.homography), scene["K"]),
                             len(inlier_indices(stats)), elapsed))
            results[name][key] = runs
    if args.output:
        settings = dict(matches=Path(args.matches).name, threshold=args.threshold, backend=args.backend,
                        min_iterations=args.min_iterations, seeds=args.seeds, pairs=keys, **extra)
        with open(args.output, "w") as f:
            json.dump(dict(settings=settings, results=results), f)
    print_results(results, keys, pairs, args.per_pair)


def print_results(results, keys, pairs, per_pair):
    paper = {s: (label, err, t) for s, label, err, t in TABLE3}
    print(f"\n{len(keys)} pairs: " + ", ".join(f"{k} ({pairs[k]['x'].shape[1]} matches)" for k in keys))
    print(f"{'Solver':<36s} {'Error':>7s} {'Error (px)':>10s} {'Inliers':>8s} {'Time (ms)':>10s}   "
          f"{'Paper: error':>12s} {'time (ms)':>9s}")
    for name, per in results.items():
        e, e_px, n, t = (np.array([run[i] for k in keys for run in per[k]]) for i in range(4))
        _, p_err, p_t = paper[name]
        print(f"{name:<36s} {e.mean():7.3f} {e_px.mean():10.3f} {n.mean():8.1f} {t.mean() * 1e3:10.1f}   "
              f"{p_err:12.3f} {p_t:9d}")
    if per_pair:
        print("\nError per pair (median over the seeds), in calibrated coordinates / in pixels")
        print(f"{'Solver':<36s} " + " ".join(f"{k:>13s}" for k in keys))
        for name, per in results.items():
            print(f"{name:<36s} " + " ".join(f"{np.median([r[0] for r in per[k]]):6.3f}/"
                                             f"{np.median([r[1] for r in per[k]]):6.3f}" for k in keys))


def compare(path_a, path_b):
    """The differences between two runs, b - a, for the pairs and solvers of both: the errors,
    averaged over the seeds of each pair, with their standard errors over the pairs; and the ratio
    of the mean times."""
    runs = []
    for path in (path_a, path_b):
        with open(path) as f:
            runs.append(json.load(f))
    (a, b) = (r["results"] for r in runs)
    warn_settings(runs[0].get("settings"), runs[1].get("settings"))
    print(f"{'Solver':<36s} {'Error':>17s} {'Error (px)':>19s} {'Time':>7s}")
    for name in sorted(set(a) & set(b)):
        keys = [k for k in a[name] if k in b[name]]
        mean = lambda res, i: np.array([np.mean([run[i] for run in res[name][k]]) for k in keys])  # noqa: E731
        d, d_px = mean(b, 0) - mean(a, 0), mean(b, 1) - mean(a, 1)
        ratio = mean(b, 3).mean() / mean(a, 3).mean()
        print(f"{name:<36s} {d.mean():+8.4f} ±{standard_error(d):6.4f} {d_px.mean():+10.4f} "
              f"±{standard_error(d_px):6.4f} {ratio:6.2f}x")


def standard_error(d):
    """The standard error of the mean of d: nan with fewer than 2 values"""
    return d.std(ddof=1) / np.sqrt(len(d)) if len(d) > 1 else np.nan


if __name__ == "__main__":
    main()
