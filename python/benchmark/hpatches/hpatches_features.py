"""
Matches of the HPatches pairs, for hpatches_benchmark.py.

For each pair (image 1 with each of the 5 others in each sequence, 580 pairs), the script matches
features and saves the matched points with the ground-truth homography. The points are centered
on each image's center (HomLib's distortion center), and the homography is in these coordinates.

Pixel coordinates: OpenCV and LightGlue put (0, 0) at the center of the top-left pixel, so the
image center is ((w - 1) / 2, (h - 1) / 2). LoMa puts (0, 0) at the top-left corner of the image
(the center of the top-left pixel is (0.5, 0.5)): the feature class shifts its points by -0.5 px.

The matches files also store the sizes of the two images of each pair. The second image's is
new: ``hpatches_benchmark.py --protocol gluefactory`` needs it, and reads it from HPatches for the
older files, which lack it.

The features and the matcher (implemented in ``benchmark/features.py``):
- ``--features sift`` (default): OpenCV SIFT, matched with Lowe's ratio test. Easy to reproduce,
  with no learned model. With OpenCV 4.8 or later, it uses SIFT's precise upscaling, which keypoint
  accuracy needs.
- ``--features superpoint``, ``disk``, ``aliked`` or ``sift`` with ``--matcher lightglue``: these
  features matched with LightGlue (Lindenberger et al., ICCV 2023), with its own extractors.
  Wadenbäck et al. (3DV 2026) evaluate on HPatches with SuperPoint and LightGlue, 1024 keypoints,
  and the shorter side of the images resized to 480 pixels: ``--features superpoint --matcher
  lightglue --num_features 1024 --resize 480`` follows their setup. The points are given in the
  original images, so the thresholds and the corner errors of hpatches_benchmark.py are in pixels of
  the original images, while glue-factory's evaluation, and probably Wadenbäck et al., measure them
  in the resized images: ``hpatches_benchmark.py --protocol gluefactory`` measures them as it does.
  LightGlue's SIFT is OpenCV SIFT without precise upscaling, which shifts its points by about
  0.1 px: the SIFTLightGlue class gives it the same SIFT with precise upscaling.
- ``--features loma``: LoMa (Nordström et al., ECCV 2026), which detects, describes and matches the
  points itself. Its pixel origin is the top-left corner of the image: the feature class shifts its
  points by -0.5 px, to OpenCV's convention.

Each match also gets a score, higher for a better match, as the ``scores`` of the robust estimators
(not sorted): 1 - the distance ratio of the ratio test for SIFT, the match confidence of LightGlue
and of LoMa.

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
import sys
from pathlib import Path

import cv2
import numpy as np
from tqdm import tqdm

from hpatches_data import download_hpatches, load_hpatches_pairs, load_matches, save_matches

# features.py is in the parent directory (benchmark/)
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from features import add_feature_args, create_matcher  # noqa: E402

SUBSET = ["v_graffiti/2", "v_boat/3", "v_wall/4", "v_bird/2",
          "i_dome/3", "i_castle/2", "i_parking/4", "i_greentea/5"]


def compute_matches(seq_dir: Path, matcher):
    pairs = {}
    for seq_name, i, path1, path2, H in tqdm(list(load_hpatches_pairs(seq_dir)), desc="Matches"):
        result = matcher(path1, path2)
        p1, p2, scores = result.points1, result.points2, result.scores
        (h1, w1), (h2, w2) = (cv2.imread(p, cv2.IMREAD_GRAYSCALE).shape for p in (path1, path2))
        c1 = np.array([(w1 - 1) / 2, (h1 - 1) / 2])
        c2 = np.array([(w2 - 1) / 2, (h2 - 1) / 2])
        T1 = np.array([[1, 0, -c1[0]], [0, 1, -c1[1]], [0, 0, 1.0]])
        T2 = np.array([[1, 0, -c2[0]], [0, 1, -c2[1]], [0, 0, 1.0]])
        H_centered = T2 @ H @ np.linalg.inv(T1)
        pairs[f"{seq_name}/{i}"] = dict(x=(p1 - c1).T, y=(p2 - c2).T, H=H_centered, scores=scores,
                                        width=w1, height=h1, width2=w2, height2=h2)
    return pairs


def main():
    parser = argparse.ArgumentParser(description="Matches of the HPatches pairs")
    parser.add_argument("--data_dir", type=str, default="data", help="Directory for HPatches data")
    add_feature_args(parser)
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
        try:
            matcher = create_matcher(args)
        except ValueError as e:
            parser.error(str(e))
        pairs = compute_matches(download_hpatches(Path(args.data_dir)), matcher)
        name = args.features if args.features == "loma" or args.matcher == "ratio" else f"{args.features}_lightglue"
        output = Path(args.output or Path(args.data_dir) / f"hpatches_{name}.npz")
        save_matches(output, pairs)
        print(f"Wrote the matches of {len(pairs)} pairs to {output}")
    if args.export_subset:
        save_matches(Path(args.export_subset), {k: pairs[k] for k in SUBSET})
        print(f"Wrote {len(SUBSET)} pairs to {args.export_subset}")


if __name__ == "__main__":
    main()
