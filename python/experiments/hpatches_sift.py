"""
HPatches benchmark of HomLib's point solvers, with OpenCV SIFT matches.

For each HPatches pair (image 1 with each of the 5 others, 580 pairs), the script matches SIFT
features (ratio test), runs LOMSAC with every point solver, and reports, for the viewpoint and
the illumination sequences:
- the share of the true inliers that LOMSAC finds (the matches within the inlier threshold of the
  ground-truth homography);
- the AUC of the corner error, mapping the image corners with the estimated model (distortion
  included) and with the ground-truth homography;
- the median time per pair.

The matches are cached in the data directory, so later runs skip the feature extraction.
``--export_subset`` writes the matches of a few pairs to a file, for the regression test
``python/tests/test_hpatches.py``.

Requirements
------------
pip install opencv-python tqdm

Usage
-----
    python hpatches_sift.py [--data_dir data] [--solvers Wadenback3DV2026TwoSided DLT] [--seeds 3]
"""

import argparse
import time
from collections import defaultdict
from pathlib import Path

import cv2
import numpy as np
from tqdm import tqdm

import homlib

HPATCHES_URL = "http://icvl.ee.ic.ac.uk/vgg/hpatches/sequences/hpatches-sequences-release.tar.gz"
HPATCHES_DIRNAME = "hpatches-sequences-release"


def download_hpatches(data_dir: Path) -> Path:
    """Download and extract HPatches if not already present."""
    seq_dir = data_dir / HPATCHES_DIRNAME
    if seq_dir.exists():
        return seq_dir
    import tarfile
    import urllib.request
    archive = data_dir / "hpatches-sequences-release.tar.gz"
    if not archive.exists():
        print(f"Downloading HPatches to {archive} ...")
        data_dir.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(HPATCHES_URL, archive)
    print("Extracting ...")
    with tarfile.open(archive, "r:gz") as tar:
        tar.extractall(data_dir)
    return seq_dir


def sift_matches(seq_dir: Path, num_features: int, ratio: float):
    """The SIFT matches of every pair, as {"seq/i": dict(x, y, H, width, height)}. The points
    are centered on each image's center, and H is the ground truth in these coordinates."""
    sift = cv2.SIFT_create(num_features)
    matcher = cv2.BFMatcher()
    pairs = {}
    for seq in tqdm(sorted(p for p in seq_dir.iterdir() if p.is_dir()), desc="SIFT"):
        im1 = cv2.imread(str(seq / "1.ppm"), cv2.IMREAD_GRAYSCALE)
        k1, d1 = sift.detectAndCompute(im1, None)
        c1 = np.array([(im1.shape[1] - 1) / 2, (im1.shape[0] - 1) / 2])
        for i in range(2, 7):
            im2 = cv2.imread(str(seq / f"{i}.ppm"), cv2.IMREAD_GRAYSCALE)
            k2, d2 = sift.detectAndCompute(im2, None)
            c2 = np.array([(im2.shape[1] - 1) / 2, (im2.shape[0] - 1) / 2])
            good = [m for m, n in matcher.knnMatch(d1, d2, k=2) if m.distance < ratio * n.distance]
            x = np.array([k1[m.queryIdx].pt for m in good]).reshape(-1, 2) - c1
            y = np.array([k2[m.trainIdx].pt for m in good]).reshape(-1, 2) - c2
            T1 = np.array([[1, 0, -c1[0]], [0, 1, -c1[1]], [0, 0, 1.0]])
            T2 = np.array([[1, 0, -c2[0]], [0, 1, -c2[1]], [0, 0, 1.0]])
            H = T2 @ np.loadtxt(seq / f"H_1_{i}") @ np.linalg.inv(T1)
            pairs[f"{seq.name}/{i}"] = dict(x=x.T, y=y.T, H=H, width=im1.shape[1], height=im1.shape[0])
    return pairs


def load_or_compute_matches(data_dir: Path, num_features: int, ratio: float):
    cache = data_dir / f"hpatches_sift_{num_features}_{ratio}.npz"
    if cache.exists():
        flat = np.load(cache)
        keys = sorted({k.rsplit(":", 1)[0] for k in flat.files})
        return {k: dict(x=flat[k + ":x"], y=flat[k + ":y"], H=flat[k + ":H"],
                        width=int(flat[k + ":size"][0]), height=int(flat[k + ":size"][1])) for k in keys}
    pairs = sift_matches(download_hpatches(data_dir), num_features, ratio)
    save_pairs(cache, pairs)
    return pairs


def save_pairs(path: Path, pairs):
    flat = {}
    for k, p in pairs.items():
        flat[k + ":x"], flat[k + ":y"], flat[k + ":H"] = p["x"], p["y"], p["H"]
        flat[k + ":size"] = np.array([p["width"], p["height"]])
    np.savez_compressed(path, **flat)


def true_inliers(pair, threshold):
    """The matches within the threshold of the ground-truth homography"""
    x, y, H = pair["x"], pair["y"], pair["H"]
    p = H @ np.vstack([x, np.ones(x.shape[1])])
    return np.flatnonzero(np.linalg.norm(p[:2] / p[2] - y, axis=0) < threshold)


def map_points(pose, p):
    """The points p of image 1 mapped to image 2 with an estimated model: undistorted with k1,
    mapped with H, distorted with k2 (division model)"""
    k1, k2 = pose.distortion_parameter, pose.distortion_parameter2
    q = np.array(pose.homography) @ np.vstack([p, 1 + k1 * np.sum(p**2, 0)])
    u = q[:2] / q[2]
    disc = 1 - 4 * k2 * np.sum(u**2, 0)
    with np.errstate(invalid="ignore"):
        return np.where(disc >= 0, u * 2 / (1 + np.sqrt(np.maximum(disc, 0))), np.nan)


def corner_error(pose, pair):
    """Mean distance between the image corners mapped by the estimated model and by the ground
    truth, in pixels"""
    w, h = pair["width"] / 2, pair["height"] / 2
    corners = np.array([[-w, w, w, -w], [-h, -h, h, h]])
    g = pair["H"] @ np.vstack([corners, np.ones(4)])
    err = np.mean(np.linalg.norm(map_points(pose, corners) - g[:2] / g[2], axis=0))
    return err if np.isfinite(err) else np.inf


def point_solvers():
    return [s for s in homlib.SolverId.__members__.values()
            if homlib.Solver(s).feature_type == homlib.FeatureType.Points]


def main():
    parser = argparse.ArgumentParser(description="HPatches + OpenCV SIFT + HomLib benchmark")
    parser.add_argument("--data_dir", type=str, default="data", help="Directory for HPatches data")
    parser.add_argument("--num_features", type=int, default=4000, help="SIFT features per image")
    parser.add_argument("--ratio", type=float, default=0.8, help="Ratio test")
    parser.add_argument("--threshold", type=float, default=3.0, help="Inlier threshold in pixels")
    parser.add_argument("--seeds", type=int, default=3, help="LOMSAC runs per pair (random seeds)")
    parser.add_argument("--solvers", nargs="*", default=None,
                        help="SolverId names (default: every point solver)")
    parser.add_argument("--export_subset", type=str, default=None,
                        help="Write the matches of a few pairs to this file, and exit")
    args = parser.parse_args()

    pairs = load_or_compute_matches(Path(args.data_dir), args.num_features, args.ratio)
    if args.export_subset:
        keys = ["v_graffiti/2", "v_boat/3", "v_wall/4", "v_bird/2",
                "i_dome/3", "i_castle/2", "i_parking/4", "i_greentea/5"]
        save_pairs(Path(args.export_subset), {k: pairs[k] for k in keys})
        print(f"Wrote {len(keys)} pairs to {args.export_subset}")
        return

    solvers = ([homlib.SolverId.__members__[s] for s in args.solvers] if args.solvers
               else point_solvers())
    results = defaultdict(lambda: defaultdict(list))
    for key, pair in tqdm(pairs.items(), desc="Pairs"):
        if pair["x"].shape[1] < 8:
            continue
        x, y = np.ascontiguousarray(pair["x"]), np.ascontiguousarray(pair["y"])
        gt = set(true_inliers(pair, args.threshold).tolist())
        kind = "viewpoint" if key.startswith("v_") else "illumination"
        for s in solvers:
            for seed in range(1, args.seeds + 1):
                lomsac = homlib.LOMSAC(homlib.Solver(s), inlier_threshold=args.threshold, random_seed=seed)
                t0 = time.perf_counter()
                pose, stats = lomsac.estimate(x, y)
                elapsed = time.perf_counter() - t0
                r = results[(s.name, kind)]
                r["recall"].append(len(gt & set(stats.inlier_indices)) / max(len(gt), 1))
                r["error"].append(corner_error(pose, pair))
                r["time"].append(elapsed)

    print(f"\n{'Solver':<46s} {'Pairs':<13s} {'True inliers':>12s}  AUC@1  AUC@3  AUC@5 AUC@10  {'Time (ms)':>9s}")
    for (name, kind), r in sorted(results.items()):
        err = np.array(r["error"])
        aucs = [np.mean(np.minimum(err, t) / t) * 100 for t in (1, 3, 5, 10)]
        aucs = [100 - a for a in aucs]   # AUC of the cumulative error curve up to t
        print(f"{name:<46s} {kind:<13s} {np.mean(r['recall']) * 100:11.2f}%"
              + "".join(f" {a:6.2f}" for a in aucs) + f"  {np.median(r['time']) * 1e3:9.2f}")


if __name__ == "__main__":
    main()
