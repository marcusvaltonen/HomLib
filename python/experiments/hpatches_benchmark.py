"""
HPatches benchmark of HomLib's point solvers, on precomputed matches (hpatches_features.py).

For each pair, the script runs a robust estimator (LOMSAC by default) with every point solver, and
reports, for the viewpoint and the illumination sequences:
- the share of the true inliers that the estimator finds, the true inliers being the matches within
  the inlier threshold of the ground-truth homography. This share ignores false inliers: the corner
  error accounts for them;
- the AUC of the corner error, which maps the image corners with the estimated model (distortion
  included) and with the ground-truth homography. Wadenbäck et al. (3DV 2026) report the AUC at
  1, 3 and 5 px on HPatches, as "mAA (AUC)";
- the median time of a run.

By default (``--protocol original``), the thresholds and the errors are in pixels of the original
images, on the 580 pairs, and the pairs with no true inlier are left out, since no solver can find
any.

``--protocol gluefactory`` follows glue-factory's HPatches evaluation (gluefactory/eval/hpatches.py
and gluefactory/datasets/hpatches.py in https://github.com/cvg/glue-factory), so that the results
compare with those measured with it:
- each image is resized to 480 px on its shorter side, the other side rounded down (as glue-factory's
  ImagePreprocessor): the points, the ground truth and the image corners are scaled to the resized
  images, so that the inlier threshold, the true inliers and the corner error are in their pixels;
- the ground truth follows glue-factory's convention, which applies the homographies of HPatches to
  coordinates with their origin at the top-left corner of the image, where hpatches_features.py
  applies them with the origin at the center of the top-left pixel (see resize_pair);
- the 8 sequences of large images that glue-factory leaves out (GLUEFACTORY_IGNORED_SCENES) are left
  out: 540 pairs;
- every pair counts in the AUC: a pair with no true inlier counts in the AUC but not in the share of
  the true inliers found, and an estimate with no inlier is a failure, with an infinite corner error;
- the table has a row for all the pairs, as glue-factory's AUCs, besides the rows for the
  illumination and the viewpoint pairs.
glue-factory integrates the AUC with the trapezoidal rule, which adds less than 50 / n points to it
(n: pairs times seeds), so less than 0.1 points on 540 pairs. The robust estimator stays HomLib's,
whereas glue-factory runs OpenCV's findHomography, with a 0.5 px threshold in its configurations of
LightGlue (``--threshold 0.5``). Older matches files, such as the regression test's subset, lack the
sizes of the second images: ``--data_dir`` then gives the directory of HPatches, whose images have
them.

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
                                 [--protocol gluefactory [--data_dir data]]
    python hpatches_benchmark.py --compare before.json after.json
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np
from tqdm import tqdm

import homlib
from hpatches_data import HPATCHES_DIRNAME, add_image_sizes, load_matches

# glue-factory's HPatches evaluation: the shorter side of the resized images, and the sequences of
# large images that it leaves out (HPatches.ignored_scenes in gluefactory/datasets/hpatches.py)
GLUEFACTORY_RESIZE = 480
GLUEFACTORY_IGNORED_SCENES = ("i_contruction", "i_crownnight", "i_dc", "i_pencils", "i_whitebuilding",
                              "v_artisans", "v_astronautis", "v_talent")


def true_inliers(pair, threshold):
    """The matches within the threshold of the ground-truth homography."""
    x, y, H = pair["x"], pair["y"], pair["H"]
    p = H @ np.vstack([x, np.ones(x.shape[1])])
    return np.flatnonzero(np.linalg.norm(p[:2] / p[2] - y, axis=0) < threshold)


def gluefactory_size(width, height, side=GLUEFACTORY_RESIZE):
    """(width, height) of an image resized as in glue-factory's evaluation (ImagePreprocessor with
    resize=side, side="short"): the shorter side becomes side, and the other is rounded down."""
    aspect_ratio = width / height
    return (int(side * aspect_ratio), side) if aspect_ratio >= 1.0 else (side, int(side / aspect_ratio))


def resize_pair(pair):
    """The pair in pixels of its two images resized as in glue-factory's evaluation. Resizing an
    image of width w to w' maps x to (x + 0.5) w' / w - 0.5 (pixel centers at integer coordinates),
    and the image center (w - 1) / 2 to (w' - 1) / 2: the points, centered on the image center,
    scale by w' / w, and likewise in y.

    The ground truth follows glue-factory's convention: glue-factory applies the homographies of
    HPatches to coordinates with their origin at the top-left corner of the image, whereas
    hpatches_features.py applies them to coordinates with their origin at the center of the
    top-left pixel, as OpenCV. Which convention HPatches itself uses is undetermined. In the
    centered coordinates, glue-factory's ground truth is T(-1/2) H T(1/2), where H is the pair's
    and T(a) the translation by (a, a). The ground truths of the two conventions map the image
    corners 0.08 px apart (median over the 540 pairs, in the resized images), but more than 1 px
    apart on 1.7% of the pairs (5 px on v_charing/3)."""
    w1, h1 = gluefactory_size(pair["width"], pair["height"])
    w2, h2 = gluefactory_size(pair["width2"], pair["height2"])
    S1 = np.diag([w1 / pair["width"], h1 / pair["height"], 1.0])
    S2 = np.diag([w2 / pair["width2"], h2 / pair["height2"], 1.0])
    H = translation(-0.5) @ pair["H"] @ translation(0.5)
    return dict(x=S1[:2, :2] @ pair["x"], y=S2[:2, :2] @ pair["y"], H=S2 @ H @ np.linalg.inv(S1),
                width=w1, height=h1, width2=w2, height2=h2)


def translation(a):
    """The homography of the translation by (a, a)"""
    return np.array([[1.0, 0, a], [0, 1, a], [0, 0, 1]])


def gluefactory_pairs(pairs, data_dir):
    """The pairs of glue-factory's evaluation, resized as in it. data_dir: the directory of HPatches
    (with hpatches-sequences-release), for the sizes of the second images, if the matches file
    lacks them; None if it has them."""
    pairs = {k: p for k, p in pairs.items() if k.split("/")[0] not in GLUEFACTORY_IGNORED_SCENES}
    if any(p["width2"] is None for p in pairs.values()):
        if data_dir is None:
            raise ValueError("the matches file lacks the sizes of the second images, which --protocol "
                             "gluefactory needs: give the directory of HPatches with --data_dir")
        add_image_sizes(pairs, Path(data_dir) / HPATCHES_DIRNAME)
    return {k: resize_pair(p) for k, p in pairs.items()}


def map_points(pose, p):
    """The points p of image 1 mapped to image 2 with an estimated model: undistorted with k1,
    mapped with H, distorted with k2 (division model); and whether each point has an image. It
    hasn't if H maps it to infinity, or if k2 has no distorted point for it (1 - 4 k2 |u|^2 < 0).
    Both are tested before the division and the square root, so that no NaN or infinity is
    created."""
    k1, k2 = pose.distortion_parameter, pose.distortion_parameter2
    # Undistorted with k1 in homogeneous coordinates (w = 1 + k1 |p|^2 as third coordinate): no
    # division by w
    q = np.array(pose.homography) @ np.vstack([p, 1 + k1 * np.sum(p**2, 0)])
    ok = q[2] != 0
    u = np.divide(q[:2], q[2], out=np.zeros_like(q[:2]), where=ok)
    disc = 1 - 4 * k2 * np.sum(u**2, 0)
    ok &= disc >= 0
    return u * 2 / (1 + np.sqrt(np.where(ok, disc, 0))), ok


def corner_error(pose, pair):
    """The mean distance between the image corners mapped by the estimated model and by the
    ground truth, in pixels; infinite if the model doesn't map every corner."""
    # The points are centered on the image center, ((w - 1) / 2, (h - 1) / 2) with OpenCV's pixel
    # centers at integer coordinates: the outer corners of the image are at (±w / 2, ±h / 2).
    w, h = pair["width"] / 2, pair["height"] / 2
    corners = np.array([[-w, w, w, -w], [-h, -h, h, h]])
    g = pair["H"] @ np.vstack([corners, np.ones(4)])
    mapped, ok = map_points(pose, corners)
    if not ok.all():
        return np.inf
    err = np.mean(np.linalg.norm(mapped - g[:2] / g[2], axis=0))
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


def main():
    parser = argparse.ArgumentParser(description="HPatches benchmark of HomLib's point solvers")
    parser.add_argument("matches", nargs="?", help="The matches file (hpatches_features.py)")
    parser.add_argument("--threshold", type=float, default=3.0,
                        help="Inlier threshold in pixels (of the resized images with --protocol gluefactory)")
    parser.add_argument("--protocol", choices=["original", "gluefactory"], default="original",
                        help="original (default): in pixels of the original images, on the 580 pairs; "
                             "gluefactory: glue-factory's HPatches evaluation, in images resized to 480 px "
                             "(shorter side), on 540 pairs")
    parser.add_argument("--data_dir", type=str, default=None,
                        help="With --protocol gluefactory, for matches files that lack the sizes of the second "
                             "images: the directory of HPatches (with hpatches-sequences-release)")
    parser.add_argument("--backend", choices=sorted(BACKENDS), default="lomsac",
                        help="Robust estimator (default: lomsac). The others need HomLib built with them")
    parser.add_argument("--seeds", type=int, default=3, help="Runs per pair (random seeds, for LOMSAC)")
    parser.add_argument("--iterations", type=int, default=None,
                        help="Run the estimator for exactly this many iterations")
    point_solvers = homlib.solver_ids(homlib.FeatureType.Points)
    parser.add_argument("--solvers", nargs="+", choices=[s.name for s in point_solvers], default=None,
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
    if args.protocol == "gluefactory":
        try:
            pairs = gluefactory_pairs(pairs, args.data_dir)
        except ValueError as e:
            parser.error(str(e))

    solvers = ([homlib.SolverId.__members__[s] for s in args.solvers] if args.solvers
               else point_solvers)
    # results[solver][pair] = per seed: (share of the true inliers found, corner error, time)
    results = {s.name: {} for s in solvers}
    for key, pair in tqdm(pairs.items(), desc="Pairs"):
        gt = set(true_inliers(pair, args.threshold).tolist())
        if not gt and args.protocol == "original":
            continue
        x, y = np.ascontiguousarray(pair["x"]), np.ascontiguousarray(pair["y"])
        for s in solvers:
            runs = []
            for seed in range(1, args.seeds + 1):
                estimator = make_estimator(args.backend, homlib.Solver(s), args.threshold, args.iterations, seed)
                t0 = time.perf_counter()
                pose, stats = estimator.estimate(x, y)
                elapsed = time.perf_counter() - t0
                found = set(inlier_indices(stats))
                # glue-factory counts a failed estimate (no inlier) as an infinite error
                failed = not found and args.protocol == "gluefactory"
                runs.append((len(gt & found) / len(gt) if gt else np.nan,
                             np.inf if failed else corner_error(pose, pair), elapsed))
            results[s.name][key] = runs
    if args.output:
        settings = dict(matches=Path(args.matches).name, threshold=args.threshold, backend=args.backend,
                        seeds=args.seeds, iterations=args.iterations, protocol=args.protocol)
        with open(args.output, "w") as f:
            # Strict JSON: null for the share of the true inliers found of a pair with none (NaN), and
            # for the corner error of a failure (infinite); load_results maps them back
            strict = {name: {key: [[v if np.isfinite(v) else None for v in run] for run in runs]
                             for key, runs in per_pair.items()} for name, per_pair in results.items()}
            json.dump(dict(settings=settings, results=strict), f, allow_nan=False)
    print_results(results, overall=args.protocol == "gluefactory")


def load_results(path):
    """(settings, results) of a file saved with --output. The files saved before the settings were
    the results alone: their settings are None. The files saved before --protocol used the original
    protocol. A null share of the true inliers found is NaN, and a null corner error is infinite
    (older files have NaN and Infinity instead)."""
    with open(path) as f:
        data = json.load(f)
    settings, results = ((None, data) if "settings" not in data
                         else ({"protocol": "original", **data["settings"]}, data["results"]))
    for per_pair in results.values():
        for runs in per_pair.values():
            for run in runs:
                run[0] = np.nan if run[0] is None else run[0]
                run[1] = np.inf if run[1] is None else run[1]
    return settings, results


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


def print_results(results, overall=False):
    """The table of the results, for the illumination and the viewpoint pairs, and with overall, for
    all the pairs (as glue-factory's evaluation)."""
    keys = next(iter(results.values()), {})
    print(f"\n{sum(kind(k) == 'viewpoint' for k in keys)} viewpoint and "
          f"{sum(kind(k) == 'illumination' for k in keys)} illumination pairs")
    print(f"{'Solver':<46s} {'Sequences':<13s} {'Inliers found':>13s}  AUC@1  AUC@3  AUC@5 AUC@10  "
          f"{'Time (ms)':>9s}")
    for name, per_pair in sorted(results.items()):
        for k in ("illumination", "viewpoint") + (("all",) if overall else ()):
            runs = [r for key, rs in per_pair.items() if k in ("all", kind(key)) for r in rs]
            if not runs:
                continue
            recall, err, t = (np.array(v, dtype=float) for v in zip(*runs))
            # The pairs with no true inlier (--protocol gluefactory) have no share of them found
            recall = recall[np.isfinite(recall)]
            # AUC of the cumulative error curve up to each threshold
            aucs = [np.mean(np.maximum(0, 1 - err / th)) * 100 for th in (1, 3, 5, 10)]
            recall = f"{np.mean(recall) * 100:12.2f}%" if len(recall) else f"{'n/a':>13s}"
            print(f"{name:<46s} {k:<13s} {recall}"
                  + "".join(f" {a:6.2f}" for a in aucs) + f"  {np.median(t) * 1e3:9.2f}")


def compare(path_a, path_b):
    """The differences between two runs, b - a, for the pairs and solvers of both: share of the
    true inliers found and AUC@5, averaged over the seeds of each pair, with their standard errors
    over the pairs (the share of the true inliers found, over the pairs with true inliers in both);
    and the ratio of the median times, which has no standard error."""
    settings_a, a = load_results(path_a)
    settings_b, b = load_results(path_b)
    warn_settings(settings_a, settings_b)
    # With two runs of glue-factory's protocol, also all the pairs together, as its AUCs
    overall = all(s is not None and s["protocol"] == "gluefactory" for s in (settings_a, settings_b))
    print(f"{'Solver':<46s} {'Sequences':<13s} {'Inliers found':>18s} {'AUC@5':>18s} {'Time':>7s}")
    for name in sorted(set(a) & set(b)):
        for k in ("illumination", "viewpoint") + (("all",) if overall else ()):
            keys = [key for key in a[name] if key in b[name] and k in ("all", kind(key))]
            if not keys:
                continue
            d_recall = (per_pair(b, name, keys, 0) - per_pair(a, name, keys, 0)) * 100
            d_recall = d_recall[np.isfinite(d_recall)]   # without the pairs with no true inlier
            d_auc = (per_pair(b, name, keys, 1) - per_pair(a, name, keys, 1)) * 100
            ratio = np.median(per_pair(b, name, keys, 2)) / np.median(per_pair(a, name, keys, 2))
            mean_recall = d_recall.mean() if len(d_recall) else np.nan
            print(f"{name:<46s} {k:<13s} {mean_recall:+8.2f} ±{standard_error(d_recall):5.2f} pts"
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
