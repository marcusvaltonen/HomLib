"""Regression test on real image pairs: LOMSAC with every point solver on 8 HPatches pairs.

The OpenCV SIFT matches of these pairs, and their ground-truth homographies, come from
``python/experiments/hpatches_sift.py --export_subset``. The test looks for them at
``$HOMLIB_TEST_DATA/hpatches_sift_subset.npz``, then downloads them from ``SUBSET_URL`` to a
cache directory. Without them, it is skipped.
"""

import hashlib
import os
import unittest
import urllib.request
from pathlib import Path

import numpy as np

import homlib

SUBSET_NAME = "hpatches_sift_subset.npz"
SUBSET_URL = os.environ.get("HOMLIB_TEST_DATA_URL", "")   # TODO: the release asset, once published
SUBSET_SHA256 = ""   # TODO: with the URL

THRESHOLD = 3.0   # inlier threshold in pixels
SEEDS = (1, 2, 3)

# The share of the true inliers that LOMSAC must find, on average over the pairs and the seeds.
# The two-sided distortion solvers are weaker on some illumination pairs, whose homography is the
# identity: there, any k1 = k2 fits the matches.
MIN_RECALL = 0.95
MIN_RECALL_TWO_SIDED = 0.7


def find_subset():
    data_dir = os.environ.get("HOMLIB_TEST_DATA")
    if data_dir and (Path(data_dir) / SUBSET_NAME).exists():
        return Path(data_dir) / SUBSET_NAME
    if not SUBSET_URL:
        return None
    cache = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "homlib" / SUBSET_NAME
    if not cache.exists():
        try:
            cache.parent.mkdir(parents=True, exist_ok=True)
            urllib.request.urlretrieve(SUBSET_URL, cache)
        except OSError:
            return None
    if SUBSET_SHA256 and hashlib.sha256(cache.read_bytes()).hexdigest() != SUBSET_SHA256:
        cache.unlink()
        return None
    return cache


def load_pairs(path):
    flat = np.load(path)
    keys = sorted({k.rsplit(":", 1)[0] for k in flat.files})
    return {k: (flat[k + ":x"], flat[k + ":y"], flat[k + ":H"]) for k in keys}


def recall(solver_id, x, y, H):
    """The share of the true inliers (within THRESHOLD of H) that LOMSAC finds, averaged over SEEDS"""
    p = H @ np.vstack([x, np.ones(x.shape[1])])
    true = set(np.flatnonzero(np.linalg.norm(p[:2] / p[2] - y, axis=0) < THRESHOLD).tolist())
    found = []
    for seed in SEEDS:
        lomsac = homlib.LOMSAC(homlib.Solver(solver_id), inlier_threshold=THRESHOLD, random_seed=seed)
        _, stats = lomsac.estimate(np.ascontiguousarray(x), np.ascontiguousarray(y))
        found.append(len(true & set(stats.inlier_indices)) / len(true))
    return np.mean(found)


class HPatchesTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = find_subset()
        if path is None:
            raise unittest.SkipTest(f"{SUBSET_NAME} not found (set HOMLIB_TEST_DATA)")
        cls.pairs = load_pairs(path)

    def test_point_solvers(self):
        for solver_id in homlib.SolverId.__members__.values():
            if homlib.Solver(solver_id).feature_type != homlib.FeatureType.Points:
                continue
            two_sided = "TwoSided" in solver_id.name
            mean = np.mean([recall(solver_id, *p) for p in self.pairs.values()])
            with self.subTest(solver=solver_id.name):
                self.assertGreater(mean, MIN_RECALL_TWO_SIDED if two_sided else MIN_RECALL)


if __name__ == "__main__":
    unittest.main()
