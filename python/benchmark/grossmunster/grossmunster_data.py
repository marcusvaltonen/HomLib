"""Grossmünster Church: the image pairs, and the matches files of grossmunster_ground_truth.py."""

import json
from pathlib import Path

import numpy as np

# The six pairs of Wadenbäck et al. (3DV 2026, Section 5), images DSC_<n>.JPG. The first five are
# those given by the authors (https://github.com/marcusvaltonen/HomLib/issues/42), who dropped
# (2515, 2516) and (2516, 2517) for too few SIFT matches; they don't know the sixth, "in the
# neighbourhood". It is here (2519, 2520), the next pair of consecutive images, with 203 matches:
# - the paper's "seven images [...] for which it was possible to extract six consecutive image
#   pairs" suggests (2515, 2517) instead, which links the seven images of the five pairs once 2516 is
#   dropped. But only about 20 of its 531 SIFT matches are on the wall (within 3 px of the ground
#   truth), and every solver finds the facades in the background instead (an error above 0.4, where
#   the paper's are all about 0.04);
# - (2511, 2512), before the first pair, has 268 matches but also about 20 on the wall, and some
#   solvers fail on it.
# grossmunster_benchmark.py --pairs runs on others.
PAIRS = [(2512, 2513), (2513, 2514), (2514, 2515), (2517, 2518), (2518, 2519), (2519, 2520)]
# The other candidate for the sixth pair, which the matches files also have
OTHER_PAIRS = [(2515, 2517)]
# The pair of Fig. 8 (inliers and rectified images), in the paper's order
FIGURE_PAIR = (2518, 2517)


def image_name(n):
    return f"DSC_{n}.JPG"


def pair_key(a, b):
    return f"{a}-{b}"


def save_matches(path: Path, pairs, scene):
    """Save the matches to a .npz file. pairs maps "a-b" to dict(x, y: 2xN points of images a and b,
    centered on the principal point; H: the ground truth in these coordinates; scores: N scores of
    the matches). scene: dict of arrays (the calibration, the camera centers, the plane points...)
    and of JSON-serializable values, saved as "scene:<name>"."""
    flat = {}
    for k, p in pairs.items():
        for name in ("x", "y", "H", "scores"):
            flat[f"{k}:{name}"] = p[name]
    for name, v in scene.items():
        flat[f"scene:{name}"] = v if isinstance(v, np.ndarray) else np.array(json.dumps(v))
    np.savez_compressed(path, **flat)


def load_matches(path: Path):
    """(pairs, scene) saved by save_matches"""
    flat = np.load(path)
    pairs, scene = {}, {}
    for f in flat.files:
        k, name = f.rsplit(":", 1)
        if k == "scene":
            v = flat[f]
            scene[name] = json.loads(v.item()) if v.dtype.kind == "U" else v
        else:
            pairs.setdefault(k, {})[name] = flat[f]
    return pairs, scene
