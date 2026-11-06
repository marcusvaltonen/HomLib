"""
Matches of the HPatches pairs, for hpatches_benchmark.py.

For each pair (image 1 with each of the 5 others in each sequence, 580 pairs), the script matches
features and saves the matched points with the ground-truth homography. The points are centered
on each image's center (HomLib's distortion center), and the homography is in these coordinates.

Pixel coordinates: OpenCV and LightGlue put (0, 0) at the center of the top-left pixel, so the
image center is ((w - 1) / 2, (h - 1) / 2). LoMa puts (0, 0) at the top-left corner of the image
(the center of the top-left pixel is (0.5, 0.5)): the script shifts its points by -0.5 px.

The features and the matcher:
- ``--features sift`` (default): OpenCV SIFT, matched with Lowe's ratio test. Easy to reproduce,
  with no learned model. With OpenCV 4.8 or later, it uses SIFT's precise upscaling, which keypoint
  accuracy needs.
- ``--features superpoint``, ``disk``, ``aliked`` or ``sift`` with ``--matcher lightglue``: these
  features matched with LightGlue (Lindenberger et al., ICCV 2023), with its own extractors.
  Wadenbäck et al. (3DV 2026) evaluate on HPatches with SuperPoint and LightGlue, 1024 keypoints,
  and the shorter side of the images resized to 480 pixels: ``--features superpoint --matcher
  lightglue --num_features 1024 --resize 480``. The points are given in the original images.
- ``--features loma``: LoMa (Nordström et al., ECCV 2026), which detects, describes and matches the
  points itself, as in hpatches_loma.py.

``--export_subset`` writes the matches of 8 pairs (4 viewpoint, 4 illumination) for the regression
test ``python/tests/test_hpatches.py``.

Requirements
------------
pip install opencv-python tqdm
For ``--matcher lightglue``: LightGlue from https://github.com/cvg/LightGlue (``pip install -e .``
in a clone). SuperPoint's weights are for non-commercial use only.
For ``--features loma``: LoMa from https://github.com/davnords/LoMa (``pip install -e .`` in a
clone; the PyPI package named ``loma`` is unrelated). Both run much faster on a GPU.

Usage
-----
    python hpatches_features.py --data_dir data --output data/hpatches_sift.npz
    python hpatches_features.py --data_dir data --features superpoint --matcher lightglue \
        --num_features 1024 --resize 480 --output data/hpatches_superpoint_lightglue.npz
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
    # enable_precise_upscale (OpenCV >= 4.8) removes a shift of the keypoints of the upscaled octave
    try:
        sift = cv2.SIFT_create(num_features, enable_precise_upscale=True)
    except TypeError:
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
        # LoMa's (0, 0) is the top-left corner of the image: shift to OpenCV's convention
        return (np.asarray(kpts1, dtype=float).reshape(-1, 2) - 0.5,
                np.asarray(kpts2, dtype=float).reshape(-1, 2) - 0.5)
    return match


def lightglue_matcher(features, num_features, resize):
    import torch
    from lightglue import ALIKED, DISK, SIFT, LightGlue, SuperPoint
    from lightglue.utils import load_image, rbd
    device = "cuda" if torch.cuda.is_available() else "mps" if torch.backends.mps.is_available() else "cpu"
    extractor_class = {"superpoint": SuperPoint, "disk": DISK, "aliked": ALIKED, "sift": SIFT}[features]
    # origin="center": (0, 0) at the center of the top-left pixel, as in OpenCV (LightGlue's default
    # for now; it plans to switch to "corner")
    extractor = extractor_class(max_num_keypoints=num_features, origin="center").eval().to(device)
    if features == "sift" and extractor.conf.backend == "opencv":
        # LightGlue's OpenCV SIFT lacks precise upscaling, which shifts the keypoints of the upscaled
        # octave: use the same SIFT with it, as sift_matcher does (OpenCV >= 4.8)
        c = extractor.conf
        extractor.sift = cv2.SIFT_create(contrastThreshold=c.detection_threshold, nfeatures=c.max_num_keypoints,
                                         edgeThreshold=c.edge_threshold, nOctaveLayers=c.num_octaves,
                                         enable_precise_upscale=True)
    matcher = LightGlue(features=features).eval().to(device)
    # resize: the shorter side of the images, in pixels (None: the original images). The extractor
    # gives the keypoints in the original images.
    options = dict(resize=resize, side="short") if resize else dict(resize=None)
    cache = {}

    def extract(path):
        if path not in cache:
            cache.clear()   # image 1 of a sequence is matched with the 5 others in a row
            cache[path] = extractor.extract(load_image(path).to(device), **options)
        return cache[path]

    @torch.no_grad()
    def match(path1, path2):
        f1, f2 = extract(path1), extract(path2)
        m = rbd(matcher({"image0": f1, "image1": f2}))["matches"]
        f1, f2 = rbd(f1), rbd(f2)
        return (f1["keypoints"][m[:, 0]].cpu().numpy().astype(float),
                f2["keypoints"][m[:, 1]].cpu().numpy().astype(float))
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
    parser.add_argument("--features", choices=["sift", "superpoint", "disk", "aliked", "loma"], default="sift")
    parser.add_argument("--matcher", choices=["ratio", "lightglue"], default="ratio",
                        help="Lowe's ratio test (OpenCV SIFT only) or LightGlue")
    parser.add_argument("--num_features", type=int, default=4000, help="Keypoints per image")
    parser.add_argument("--ratio", type=float, default=None, help="Ratio test threshold (default: 0.8)")
    parser.add_argument("--resize", type=int, default=None,
                        help="With LightGlue: resize the shorter side of the images to this size")
    parser.add_argument("--output", type=str, default=None,
                        help="The matches file to write (default: <data_dir>/hpatches_<features>[_lightglue].npz)")
    parser.add_argument("--matches", type=str, default=None, help="An existing matches file (with --export_subset)")
    parser.add_argument("--export_subset", type=str, default=None,
                        help="Write the matches of the 8 pairs of the regression test to this file")
    args = parser.parse_args()

    if args.matches:
        if not args.export_subset:
            parser.error("--matches needs --export_subset")
        pairs = load_matches(Path(args.matches))
    else:
        if args.matcher == "lightglue" and args.features == "loma":
            parser.error("LoMa matches its own points: --features loma doesn't take --matcher lightglue")
        if args.resize and args.matcher != "lightglue":
            parser.error("--resize needs --matcher lightglue")
        if args.ratio is not None and args.matcher != "ratio":
            parser.error("--ratio needs --matcher ratio")
        if args.features == "loma":
            match = loma_matcher(args.num_features)
        elif args.matcher == "lightglue":
            match = lightglue_matcher(args.features, args.num_features, args.resize)
        elif args.features == "sift":
            match = sift_matcher(args.num_features, 0.8 if args.ratio is None else args.ratio)
        else:
            parser.error(f"--features {args.features} needs --matcher lightglue")
        pairs = compute_matches(download_hpatches(Path(args.data_dir)), match)
        name = args.features if args.features == "loma" or args.matcher == "ratio" else f"{args.features}_lightglue"
        output = Path(args.output or Path(args.data_dir) / f"hpatches_{name}.npz")
        save_matches(output, pairs)
        print(f"Wrote the matches of {len(pairs)} pairs to {output}")
    if args.export_subset:
        save_matches(Path(args.export_subset), {k: pairs[k] for k in SUBSET})
        print(f"Wrote {len(SUBSET)} pairs to {args.export_subset}")


if __name__ == "__main__":
    main()
