"""Tests for sampler selection and score/probability passing across all backends."""
import unittest

import numpy as np
import homlib

from tests.helpers import normalize_H, retry


def _make_scene(seed=42, n_inliers=80, n_outliers=20):
    rng = np.random.RandomState(seed)
    H = np.array([
        [1.2, 0.1, 0.05],
        [0.05, 0.9, -0.02],
        [0.001, -0.001, 1.0]
    ])
    N = n_inliers + n_outliers
    x = rng.randn(2, N) * 0.3
    y_hom = H @ np.vstack((x, np.ones((1, N))))
    y = y_hom[:2] / y_hom[2]
    x += rng.randn(2, N) * 1e-4
    y += rng.randn(2, N) * 1e-4
    x[:, :n_outliers] = rng.randn(2, n_outliers) * 10
    y[:, :n_outliers] = rng.randn(2, n_outliers) * 10

    inlier_mask = np.zeros(N, dtype=bool)
    inlier_mask[n_outliers:] = True
    scores = np.where(inlier_mask, rng.uniform(0.7, 1.0, N), rng.uniform(0.0, 0.3, N))
    order = np.argsort(-scores)
    x = np.ascontiguousarray(x[:, order])
    y = np.ascontiguousarray(y[:, order])
    scores = scores[order]
    return x, y, scores, H, n_outliers


# ---------------------------------------------------------------------------
# GC-RANSAC sampler tests
# ---------------------------------------------------------------------------

@unittest.skipUnless(homlib._HAS_GCRANSAC, "GC-RANSAC not built")
class GCRANSACSamplerTestCase(unittest.TestCase):
    def setUp(self):
        self.x, self.y, self.scores, self.H, self.n_outliers = _make_scene()
        self.N = self.x.shape[1]

    @retry()
    def test_prosac_sampler(self):
        solver = homlib.Wadenback3DV2026OneSided()
        gc = homlib.GCRANSAC(solver, threshold=0.07,
                             sampler=homlib.GCSamplerType.PROSAC)
        est, stats = gc.estimate(self.x, self.y, scores=self.scores)
        H_est = normalize_H(est.homography)
        H_gt = normalize_H(self.H)
        self.assertLess(np.linalg.norm(H_est - H_gt), 0.1)
        self.assertGreater(len(stats.inliers), self.N - self.n_outliers - 10)

    @retry()
    def test_arsampler(self):
        solver = homlib.Wadenback3DV2026OneSided()
        gc = homlib.GCRANSAC(solver, threshold=0.07,
                             sampler=homlib.GCSamplerType.ARSampler)
        est, stats = gc.estimate(self.x, self.y, scores=self.scores)
        H_est = normalize_H(est.homography)
        H_gt = normalize_H(self.H)
        self.assertLess(np.linalg.norm(H_est - H_gt), 0.1)
        self.assertGreater(len(stats.inliers), self.N - self.n_outliers - 10)

    @retry()
    def test_importance_sampler(self):
        solver = homlib.Wadenback3DV2026OneSided()
        gc = homlib.GCRANSAC(solver, threshold=0.07,
                             sampler=homlib.GCSamplerType.ImportanceSampler)
        est, stats = gc.estimate(self.x, self.y, scores=self.scores)
        H_est = normalize_H(est.homography)
        H_gt = normalize_H(self.H)
        self.assertLess(np.linalg.norm(H_est - H_gt), 0.1)

    def test_arsampler_no_scores_raises(self):
        solver = homlib.Wadenback3DV2026OneSided()
        gc = homlib.GCRANSAC(solver, threshold=0.07,
                             sampler=homlib.GCSamplerType.ARSampler)
        with self.assertRaises(RuntimeError):
            gc.estimate(self.x, self.y)

    def test_importance_sampler_no_scores_raises(self):
        solver = homlib.Wadenback3DV2026OneSided()
        gc = homlib.GCRANSAC(solver, threshold=0.07,
                             sampler=homlib.GCSamplerType.ImportanceSampler)
        with self.assertRaises(RuntimeError):
            gc.estimate(self.x, self.y)

    def test_namespace_aliases(self):
        self.assertIs(homlib.GCRANSAC.SamplerType, homlib.GCSamplerType)
        self.assertIs(homlib.GCRANSAC.Settings, homlib.GCRANSACSettings)
        self.assertIs(homlib.GCRANSAC.Statistics, homlib.GCRANSACStatistics)
        self.assertIs(homlib.GCRANSAC.NonMinimalStrategy, homlib.NonMinimalStrategy)


# ---------------------------------------------------------------------------
# MAGSAC sampler tests
# ---------------------------------------------------------------------------

@unittest.skipUnless(homlib._HAS_MAGSAC, "MAGSAC not built")
@unittest.skipUnless(homlib._HAS_GCRANSAC, "GCSamplerType requires GC-RANSAC build")
class MAGSACSamplerTestCase(unittest.TestCase):
    def setUp(self):
        self.x, self.y, self.scores, self.H, self.n_outliers = _make_scene()
        self.N = self.x.shape[1]

    @retry()
    def test_prosac_sampler(self):
        solver = homlib.Wadenback3DV2026OneSided()
        mg = homlib.MAGSAC(solver, sigma_max=0.07,
                           sampler=homlib.GCSamplerType.PROSAC)
        est, stats = mg.estimate(self.x, self.y, scores=self.scores)
        H_est = normalize_H(est.homography)
        H_gt = normalize_H(self.H)
        self.assertLess(np.linalg.norm(H_est - H_gt), 0.1)

    @retry()
    def test_arsampler(self):
        solver = homlib.Wadenback3DV2026OneSided()
        mg = homlib.MAGSAC(solver, sigma_max=0.07,
                           sampler=homlib.GCSamplerType.ARSampler)
        est, stats = mg.estimate(self.x, self.y, scores=self.scores)
        H_est = normalize_H(est.homography)
        H_gt = normalize_H(self.H)
        self.assertLess(np.linalg.norm(H_est - H_gt), 0.1)

    def test_arsampler_no_scores_raises(self):
        solver = homlib.Wadenback3DV2026OneSided()
        mg = homlib.MAGSAC(solver, sigma_max=0.07,
                           sampler=homlib.GCSamplerType.ARSampler)
        with self.assertRaises(RuntimeError):
            mg.estimate(self.x, self.y)


# ---------------------------------------------------------------------------
# SupeRANSAC sampler tests
# ---------------------------------------------------------------------------

@unittest.skipUnless(homlib._HAS_SUPERANSAC, "SupeRANSAC not built")
class SuperansacSamplerTestCase(unittest.TestCase):
    def setUp(self):
        self.x, self.y, self.scores, self.H, self.n_outliers = _make_scene()
        self.N = self.x.shape[1]

    def test_arsampler_with_scores(self):
        solver = homlib.Wadenback3DV2026OneSided()
        sa = homlib.SupeRANSAC(solver, inlier_threshold=0.07,
                                sampler=homlib.SuperansacSamplerType.ARSampler,
                                local_optimization=homlib.SuperansacLocalOptType.None_,
                                final_optimization=homlib.SuperansacLocalOptType.None_)
        est, stats = sa.estimate(self.x, self.y, scores=self.scores)
        H_est = normalize_H(est.homography)
        H_gt = normalize_H(self.H)
        self.assertLess(np.linalg.norm(H_est - H_gt), 0.1)
        self.assertGreater(len(stats.inliers), self.N - self.n_outliers - 10)

    def test_importance_sampler_with_scores(self):
        solver = homlib.Wadenback3DV2026OneSided()
        sa = homlib.SupeRANSAC(solver, inlier_threshold=0.07,
                                sampler=homlib.SuperansacSamplerType.ImportanceSampler,
                                local_optimization=homlib.SuperansacLocalOptType.None_,
                                final_optimization=homlib.SuperansacLocalOptType.None_)
        est, stats = sa.estimate(self.x, self.y, scores=self.scores)
        H_est = normalize_H(est.homography)
        H_gt = normalize_H(self.H)
        self.assertLess(np.linalg.norm(H_est - H_gt), 0.1)

    def test_arsampler_no_scores_raises(self):
        solver = homlib.Wadenback3DV2026OneSided()
        sa = homlib.SupeRANSAC(solver, inlier_threshold=0.07,
                                sampler=homlib.SuperansacSamplerType.ARSampler,
                                local_optimization=homlib.SuperansacLocalOptType.None_,
                                final_optimization=homlib.SuperansacLocalOptType.None_)
        with self.assertRaises(RuntimeError):
            sa.estimate(self.x, self.y)

    def test_namespace_aliases(self):
        self.assertIs(homlib.SupeRANSAC.SamplerType, homlib.SuperansacSamplerType)
        self.assertIs(homlib.SupeRANSAC.Settings, homlib.SuperansacSettings)
        self.assertIs(homlib.SupeRANSAC.ScoringType, homlib.SuperansacScoringType)
        self.assertIs(homlib.SupeRANSAC.LocalOptType, homlib.SuperansacLocalOptType)

    def test_all_settings_exposed(self):
        s = homlib.SuperansacSettings()
        self.assertTrue(hasattr(s, 'ar_sampler_settings'))
        self.assertTrue(hasattr(s, 'lo_settings'))
        self.assertTrue(hasattr(s, 'final_opt_settings'))
        self.assertTrue(hasattr(s, 'degensac_settings'))
        self.assertTrue(hasattr(s, 'local_optimization_inside_loop'))
        self.assertTrue(hasattr(s, 'fundamental_orientation_check'))
        self.assertTrue(hasattr(s, 'homography_bundle_refinement'))
        self.assertTrue(hasattr(s, 'final_refinement_rounds'))
        self.assertTrue(hasattr(s, 'essential_cheirality_point_number'))

    def test_ar_sampler_settings_writable(self):
        s = homlib.SuperansacSettings()
        s.ar_sampler_settings.estimator_variance = 0.5
        self.assertAlmostEqual(s.ar_sampler_settings.estimator_variance, 0.5)
        s.ar_sampler_settings.randomness = 0.02
        self.assertAlmostEqual(s.ar_sampler_settings.randomness, 0.02)

    def test_lo_settings_writable(self):
        s = homlib.SuperansacSettings()
        s.lo_settings.max_iterations = 100
        self.assertEqual(s.lo_settings.max_iterations, 100)
        s.lo_settings.spatial_coherence_weight = 0.5
        self.assertAlmostEqual(s.lo_settings.spatial_coherence_weight, 0.5)
