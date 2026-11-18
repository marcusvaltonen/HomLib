"""Tests of the benchmarks' --sampler (python/benchmark/common.py): the samplers of each backend, the
matches sorted by decreasing score for PROSAC and its inliers mapped back to their original indices,
and PROSAC on the 8 HPatches pairs of test_hpatches.py."""

import sys
import unittest
from pathlib import Path

import numpy as np

import homlib
from tests.helpers import fetch_test_data
from tests.test_hpatches import SUBSET_NAME, SUBSET_SHA256, SUBSET_URL, THRESHOLD, load_pairs

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmark"))
from common import (  # noqa: E402
    SAMPLERS, check_sampler, make_estimator, original_indices, prepare_matches, with_sampler)


class PrepareMatchesTestCase(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(0)
        self.x, self.y = rng.normal(size=(2, 7)), rng.normal(size=(2, 7))
        self.scores = np.array([0.2, 0.9, 0.5, 0.9, 0.1, 0.5, 0.7])

    def test_prosac_sorted_stable(self):
        x, y, scores, order = prepare_matches(self.x, self.y, self.scores, "prosac")
        # Decreasing scores; equal scores (0.9, 0.5) in their original order
        np.testing.assert_array_equal(order, [1, 3, 6, 2, 5, 0, 4])
        np.testing.assert_array_equal(x, self.x[:, order])
        np.testing.assert_array_equal(y, self.y[:, order])
        self.assertIsNone(scores)   # PROSAC uses only the order
        self.assertTrue(x.flags.c_contiguous and y.flags.c_contiguous)

    def test_original_indices(self):
        x, _, _, order = prepare_matches(self.x, self.y, self.scores, "prosac")
        inliers = [0, 2, 5]   # indices in the sorted matches
        np.testing.assert_array_equal(self.x[:, original_indices(inliers, order)], x[:, inliers])
        np.testing.assert_array_equal(original_indices(inliers, None), inliers)
        self.assertEqual(len(original_indices([], order)), 0)

    def test_unchanged_order(self):
        for sampler in ("uniform", "importance", "ar"):
            with self.subTest(sampler=sampler):
                x, y, scores, order = prepare_matches(self.x, self.y, self.scores, sampler)
                self.assertIsNone(order)
                np.testing.assert_array_equal(x, self.x)
                np.testing.assert_array_equal(y, self.y)
                if sampler == "uniform":
                    self.assertIsNone(scores)
                else:
                    np.testing.assert_array_equal(scores, self.scores)

    def test_no_scores(self):
        self.assertIsNone(prepare_matches(self.x, self.y, None, "uniform")[2])
        for sampler in ("prosac", "importance", "ar"):
            with self.subTest(sampler=sampler), self.assertRaises(ValueError):
                prepare_matches(self.x, self.y, None, sampler)

    def test_check_sampler(self):
        check_sampler("lomsac", "uniform")
        for backend in ("gcransac", "magsac", "superansac"):
            for sampler in ("uniform", "prosac", "importance", "ar"):
                check_sampler(backend, sampler)
        for backend, sampler in (("lomsac", "prosac"), ("lomsac", "ar"), ("gcransac", "napsac")):
            with self.subTest(backend=backend, sampler=sampler), self.assertRaises(ValueError):
                check_sampler(backend, sampler)

    def test_with_sampler(self):
        self.assertIsNone(with_sampler(None))
        self.assertEqual(with_sampler({"backend": "gcransac"})["sampler"], "uniform")
        self.assertEqual(with_sampler({"backend": "superansac"})["sampler"], "prosac (unsorted)")
        self.assertEqual(with_sampler({"backend": "superansac"}, "prosac")["sampler"], "prosac")
        self.assertEqual(with_sampler({"backend": "superansac", "sampler": "ar"})["sampler"], "ar")


BUILT = {"gcransac": homlib._HAS_GCRANSAC, "magsac": homlib._HAS_MAGSAC, "superansac": homlib._HAS_SUPERANSAC}


class MakeEstimatorTestCase(unittest.TestCase):
    def test_sampler_setting(self):
        """Uniform by default for every backend (SupeRANSAC's own default is PROSAC)"""
        for backend, built in BUILT.items():
            if not built:
                continue
            enum = homlib.SuperansacSamplerType if backend == "superansac" else homlib.GCSamplerType
            for sampler, name in SAMPLERS[backend].items():
                with self.subTest(backend=backend, sampler=sampler):
                    est = make_estimator(backend, homlib.Solver("DLT"), 3.0, None, 1, sampler=sampler)
                    self.assertEqual(est.settings.sampler, getattr(enum, name))
            with self.subTest(backend=backend, sampler="default"):
                est = make_estimator(backend, homlib.Solver("DLT"), 3.0, None, 1)
                self.assertEqual(est.settings.sampler, enum.Uniform)


class HPatchesProsacTestCase(unittest.TestCase):
    """PROSAC through prepare_matches on the 8 HPatches pairs, with scores that decrease with the
    distance to the ground truth (the subset has no scores of its own): the inliers mapped back to
    the original indices must be the true inliers."""

    @classmethod
    def setUpClass(cls):
        cls.pairs = load_pairs(fetch_test_data(SUBSET_NAME, SUBSET_URL, SUBSET_SHA256))

    def test_prosac_finds_true_inliers(self):
        backends = [b for b, built in BUILT.items() if built]
        if not backends:
            self.skipTest("no backend with PROSAC built")
        for backend in backends:
            with self.subTest(backend=backend):
                shares = []
                for x, y, H in self.pairs.values():
                    p = H @ np.vstack([x, np.ones(x.shape[1])])
                    distance = np.linalg.norm(p[:2] / p[2] - y, axis=0)
                    true = set(np.flatnonzero(distance < THRESHOLD).tolist())
                    xs, ys, scores, order = prepare_matches(x, y, 1 / (1 + distance), "prosac")
                    est = make_estimator(backend, homlib.Solver("DLT"), THRESHOLD, None, 1, sampler="prosac")
                    _, stats = est.estimate(xs, ys, scores=scores)
                    found = stats.inlier_indices if hasattr(stats, "inlier_indices") else stats.inliers
                    found = set(original_indices(found, order).tolist())
                    shares.append(len(true & found) / len(true))
                self.assertGreater(np.mean(shares), 0.95)


if __name__ == "__main__":
    unittest.main()
