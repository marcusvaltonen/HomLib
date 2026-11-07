"""Regression test on real image pairs: LOMSAC with every point solver on 8 HPatches pairs.

The OpenCV SIFT matches of these pairs, and their ground-truth homographies, come from
``python/experiments/hpatches_features.py --export_subset``. The test gets them with
``tests.helpers.fetch_test_data``: from ``$HOMLIB_TEST_DATA``, or downloaded from ``SUBSET_URL``, and
checked against ``SUBSET_SHA256``. Without them, it is skipped, except in CI
(``HOMLIB_REQUIRE_TEST_DATA``), where it fails.
"""

import os
import unittest

import numpy as np

import homlib
from tests.helpers import fetch_test_data

SUBSET_NAME = "hpatches_sift_subset.npz"
# The release asset, and the SHA-256 of the file that the test expects: a new version of the data
# changes both (see docs/source/guides/benchmarking.rst, "Test data")
SUBSET_URL = os.environ.get(
    "HOMLIB_TEST_DATA_URL",
    "https://github.com/marcusvaltonen/HomLib/releases/download/test-data-v1/" + SUBSET_NAME)
SUBSET_SHA256 = "699b605f421efedd1bba58eace3398cba645701c3a4a2492b0a949d5e6037138"

THRESHOLD = 3.0   # inlier threshold in pixels
SEEDS = (1, 2, 3)

# The share of the true inliers that each solver must find in LOMSAC, on average over the pairs and
# the seeds: a little below the measured values, with a margin for the random number generators,
# which differ between platforms. The Kukelova and Wadenbäck two-sided solvers are weaker on some
# illumination pairs, whose homography is the identity: there, any k1 = k2 fits the matches.
MIN_RECALL = 0.95
MIN_RECALL_OF = {
    "KukelovaCVPR2015TwoSided": 0.8,
    "KukelovaCVPR2015TwoSided6Pt": 0.75,
    "KukelovaCVPR2015TwoSidedEqual": 0.8,
    "KukelovaCVPR2015TwoSidedEqual6Pt": 0.75,
    "Wadenback3DV2026TwoSided": 0.85,
    "Wadenback3DV2026TwoSidedEqual": 0.85,
}


def load_pairs(path):
    flat = np.load(path)
    keys = sorted({k.rsplit(":", 1)[0] for k in flat.files})
    return {k: (flat[k + ":x"], flat[k + ":y"], flat[k + ":H"]) for k in keys}


def recall(solver_id, x, y, H):
    """The share of the true inliers (within THRESHOLD of H) that LOMSAC finds, over SEEDS."""
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
        cls.pairs = load_pairs(fetch_test_data(SUBSET_NAME, SUBSET_URL, SUBSET_SHA256))

    def test_point_solvers(self):
        for solver_id in homlib.SolverId.__members__.values():
            if homlib.Solver(solver_id).feature_type != homlib.FeatureType.Points:
                continue
            with self.subTest(solver=solver_id.name):
                mean = np.mean([recall(solver_id, *p) for p in self.pairs.values()])
                self.assertGreater(mean, MIN_RECALL_OF.get(solver_id.name, MIN_RECALL))


if __name__ == "__main__":
    unittest.main()
