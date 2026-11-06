"""
Matches of the HPatches pairs, for hpatches_benchmark.py.

For each pair (image 1 with each of the 5 others in each sequence, 580 pairs), the script matches
features and saves the matched points with the ground-truth homography. The points are centered
on each image's center (HomLib's distortion center), and the homography is in these coordinates.

Features:
- ``sift`` (default): OpenCV SIFT with Lowe's ratio test. Easy to reproduce, no learned model.
- ``loma``: LoMa (Nordström et al., ECCV 2026), as in hpatches_loma.py.

``--export_subset`` writes the matches of 8 pairs (4 viewpoint, 4 illumination) for the regression
test ``python/tests/test_hpatches.py``.

Requirements
------------
pip install opencv-python tqdm
For ``--features loma``: LoMa from https://github.com/davnords/LoMa (``pip install -e .`` in a
clone; the PyPI package named ``loma`` is unrelated).

Usage
-----
    python hpatches_features.py --data_dir data --output data/hpatches_sift.npz
    python hpatches_features.py --data_dir data --features loma --output data/hpatches_loma.npz
    python hpatches_features.py --matches data/hpatches_sift.npz --export_subset hpatches_sift_subset.npz
"""

import argparse
from pathlib import Path

import cv2
import numpy as np
from tqdm import tqdm

from hpatches_data import download_hpatches, load_hpatches_pairs, load_matches, save_matches

SUBSET = ["v_graffiti/2", "v_boat/3", "v_wall/4", "v_bird/2",
          "i_dome/3", "i_castle/2", "i_parking/4", "i_greentea/5"]


def sift_matcher(num_features, ratio):
    """match(img1_path, img2_path) -> (points1, points2), Nx2, in pixels."""
    sift = cv2.SIFT_create(num_features)
    bf = cv2.BFMatcher()
    cache = {}

    def detect(path):
        if path not in cache:
            cache.clear()   # image 1 of a sequence is matched with the 5 others in a row
            cache[path] = sift.detectAndCompute(cv2.imread(path, cv2.IMREAD_GRAYSCALE), None)
        return cache[path]

    def match(path1, path2):
        k1, d1 = detect(path1)
        k2, d2 = sift.detectAndCompute(cv2.imread(path2, cv2.IMREAD_GRAYSCALE), None)
        good = [m for m, n in bf.knnMatch(d1, d2, k=2) if m.distance < ratio * n.distance]
        return (np.array([k1[m.queryIdx].pt for m in good]).reshape(-1, 2),
                np.array([k2[m.trainIdx].pt for m in good]).reshape(-1, 2))
    return match


def loma_matcher(num_keypoints):
    from loma import LoMa, LoMaB
    model = LoMa(LoMaB(num_keypoints=num_keypoints))

    def match(path1, path2):
        kpts1, kpts2 = model.match(path1, path2)
        return np.asarray(kpts1, dtype=float).reshape(-1, 2), np.asarray(kpts2, dtype=float).reshape(-1, 2)
    return match


def compute_matches(seq_dir: Path, match):
    pairs = {}
    for seq_name, i, path1, path2, H in tqdm(list(load_hpatches_pairs(seq_dir)), desc="Matches"):
        p1, p2 = match(path1, path2)
        (h1, w1), (h2, w2) = (cv2.imread(p, cv2.IMREAD_GRAYSCALE).shape for p in (path1, path2))
        c1 = np.array([(w1 - 1) / 2, (h1 - 1) / 2])
        c2 = np.array([(w2 - 1) / 2, (h2 - 1) / 2])
        T1 = np.array([[1, 0, -c1[0]], [0, 1, -c1[1]], [0, 0, 1.0]])
        T2 = np.array([[1, 0, -c2[0]], [0, 1, -c2[1]], [0, 0, 1.0]])
        H_centered = T2 @ H @ np.linalg.inv(T1)
        pairs[f"{seq_name}/{i}"] = dict(x=(p1 - c1).T, y=(p2 - c2).T, H=H_centered, width=w1, height=h1)
    return pairs


def main():
    parser = argparse.ArgumentParser(description="Matches of the HPatches pairs")
    parser.add_argument("--data_dir", type=str, default="data", help="Directory for HPatches data")
    parser.add_argument("--features", choices=["sift", "loma"], default="sift")
    parser.add_argument("--num_features", type=int, default=4000, help="SIFT features, or LoMa keypoints")
    parser.add_argument("--ratio", type=float, default=0.8, help="Ratio test (SIFT)")
    parser.add_argument("--output", type=str, default=None,
                        help="The matches file to write (default: <data_dir>/hpatches_<features>.npz)")
    parser.add_argument("--matches", type=str, default=None, help="An existing matches file (with --export_subset)")
    parser.add_argument("--export_subset", type=str, default=None,
                        help="Write the matches of the 8 pairs of the regression test to this file")
    args = parser.parse_args()

    if args.matches:
        if not args.export_subset:
            parser.error("--matches needs --export_subset")
        pairs = load_matches(Path(args.matches))
    else:
        match = (sift_matcher(args.num_features, args.ratio) if args.features == "sift"
                 else loma_matcher(args.num_features))
        pairs = compute_matches(download_hpatches(Path(args.data_dir)), match)
        output = Path(args.output or Path(args.data_dir) / f"hpatches_{args.features}.npz")
        save_matches(output, pairs)
        print(f"Wrote the matches of {len(pairs)} pairs to {output}")
    if args.export_subset:
        save_matches(Path(args.export_subset), {k: pairs[k] for k in SUBSET})
        print(f"Wrote {len(SUBSET)} pairs to {args.export_subset}")


if __name__ == "__main__":
    main()
