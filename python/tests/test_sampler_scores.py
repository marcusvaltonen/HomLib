"""Tests for sampler selection, score passing, settings, and aliases."""
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
    scores = np.where(
        inlier_mask,
        rng.uniform(0.7, 1.0, N),
        rng.uniform(0.0, 0.3, N),
    )
    order = np.argsort(-scores)
    x = np.ascontiguousarray(x[:, order])
    y = np.ascontiguousarray(y[:, order])
    scores = scores[order]
    return x, y, scores, H, n_outliers


def _assert_good(tc, est, stats, H_gt, N, n_out,
                 check_inliers=True):
    H_est = normalize_H(est.homography)
    H_gt_n = normalize_H(H_gt)
    tc.assertLess(np.linalg.norm(H_est - H_gt_n), 0.1)
    if check_inliers:
        tc.assertGreater(len(stats.inliers), N - n_out - 10)


# ---------------------------------------------------------------
# GC-RANSAC
# ---------------------------------------------------------------

@unittest.skipUnless(homlib._HAS_GCRANSAC, "GC-RANSAC not built")
class GCRANSACSamplerTestCase(unittest.TestCase):
    def setUp(self):
        scene = _make_scene()
        self.x, self.y, self.scores, self.H, self.n_out = scene
        self.N = self.x.shape[1]
        self.solver = homlib.Wadenback3DV2026OneSided()
        self.thr = 0.07

    # --- Samplers with scores ---

    @retry()
    def test_uniform_no_scores(self):
        gc = homlib.GCRANSAC(
            self.solver, threshold=self.thr)
        est, stats = gc.estimate(self.x, self.y)
        _assert_good(
            self, est, stats, self.H, self.N, self.n_out)

    @retry()
    def test_uniform_with_scores_ignored(self):
        gc = homlib.GCRANSAC(
            self.solver, threshold=self.thr,
            sampler=homlib.GCSamplerType.Uniform)
        est, stats = gc.estimate(
            self.x, self.y, scores=self.scores)
        _assert_good(
            self, est, stats, self.H, self.N, self.n_out)

    @retry()
    def test_prosac_real_scores(self):
        gc = homlib.GCRANSAC(
            self.solver, threshold=self.thr,
            sampler=homlib.GCSamplerType.PROSAC)
        est, stats = gc.estimate(
            self.x, self.y, scores=self.scores)
        _assert_good(
            self, est, stats, self.H, self.N, self.n_out)

    @retry()
    def test_prosac_linear_scores(self):
        scores_lin = np.linspace(1.0, 0.01, self.N)
        gc = homlib.GCRANSAC(
            self.solver, threshold=self.thr,
            sampler=homlib.GCSamplerType.PROSAC)
        est, stats = gc.estimate(
            self.x, self.y, scores=scores_lin)
        _assert_good(
            self, est, stats, self.H, self.N, self.n_out)

    @retry()
    def test_prosac_no_scores(self):
        gc = homlib.GCRANSAC(
            self.solver, threshold=self.thr,
            sampler=homlib.GCSamplerType.PROSAC)
        est, stats = gc.estimate(self.x, self.y)
        _assert_good(
            self, est, stats, self.H, self.N, self.n_out)

    @retry()
    def test_arsampler_real_scores(self):
        gc = homlib.GCRANSAC(
            self.solver, threshold=self.thr,
            sampler=homlib.GCSamplerType.ARSampler)
        est, stats = gc.estimate(
            self.x, self.y, scores=self.scores)
        _assert_good(
            self, est, stats, self.H, self.N, self.n_out)

    @retry()
    def test_importance_sampler_real_scores(self):
        gc = homlib.GCRANSAC(
            self.solver, threshold=self.thr,
            sampler=homlib.GCSamplerType.ImportanceSampler)
        est, stats = gc.estimate(
            self.x, self.y, scores=self.scores)
        _assert_good(
            self, est, stats, self.H, self.N, self.n_out,
            check_inliers=False)

    # --- Error on missing scores ---

    def test_arsampler_no_scores_raises(self):
        gc = homlib.GCRANSAC(
            self.solver, threshold=self.thr,
            sampler=homlib.GCSamplerType.ARSampler)
        with self.assertRaises(RuntimeError):
            gc.estimate(self.x, self.y)

    def test_importance_sampler_no_scores_raises(self):
        gc = homlib.GCRANSAC(
            self.solver, threshold=self.thr,
            sampler=homlib.GCSamplerType.ImportanceSampler)
        with self.assertRaises(RuntimeError):
            gc.estimate(self.x, self.y)

    # --- NonMinimalStrategy ---

    @retry()
    def test_non_minimal_lomsac(self):
        nms = homlib.NonMinimalStrategy.Lomsac
        gc = homlib.GCRANSAC(
            self.solver, threshold=self.thr,
            non_minimal_strategy=nms)
        est, stats = gc.estimate(self.x, self.y)
        _assert_good(
            self, est, stats, self.H, self.N, self.n_out,
            check_inliers=False)

    @retry()
    def test_non_minimal_multisample(self):
        nms = homlib.NonMinimalStrategy.MultiSample
        gc = homlib.GCRANSAC(
            self.solver, threshold=self.thr,
            non_minimal_strategy=nms)
        est, stats = gc.estimate(self.x, self.y)
        _assert_good(
            self, est, stats, self.H, self.N, self.n_out)

    # --- Settings ---

    def test_settings_defaults(self):
        s = homlib.GCRANSACSettings()
        self.assertEqual(s.max_iterations, 5000)
        self.assertEqual(s.min_iterations, 50)
        self.assertAlmostEqual(s.confidence, 0.99)
        self.assertTrue(s.do_local_optimization)

    def test_settings_sampler_field(self):
        s = homlib.GCRANSACSettings()
        s.sampler = homlib.GCSamplerType.PROSAC
        self.assertEqual(
            s.sampler, homlib.GCSamplerType.PROSAC)

    def test_settings_non_minimal_strategy_field(self):
        s = homlib.GCRANSACSettings()
        nms = homlib.NonMinimalStrategy.Lomsac
        s.non_minimal_strategy = nms
        self.assertEqual(s.non_minimal_strategy, nms)

    # --- Namespace aliases ---

    def test_namespace_aliases(self):
        self.assertIs(
            homlib.GCRANSAC.SamplerType,
            homlib.GCSamplerType)
        self.assertIs(
            homlib.GCRANSAC.Settings,
            homlib.GCRANSACSettings)
        self.assertIs(
            homlib.GCRANSAC.Statistics,
            homlib.GCRANSACStatistics)
        self.assertIs(
            homlib.GCRANSAC.NonMinimalStrategy,
            homlib.NonMinimalStrategy)

    # --- Enum completeness ---

    def test_sampler_enum_values(self):
        members = {e.name for e in homlib.GCSamplerType}
        expected = {
            "Uniform", "PROSAC",
            "ImportanceSampler", "ARSampler",
        }
        self.assertEqual(members, expected)

    def test_non_minimal_enum_values(self):
        members = {e.name for e in homlib.NonMinimalStrategy}
        expected = {
            "Lomsac", "SelectAndRefine", "MultiSample",
        }
        self.assertEqual(members, expected)

    def test_non_minimal_enum_iterable(self):
        values = list(homlib.NonMinimalStrategy)
        self.assertEqual(len(values), 3)
        self.assertIn(homlib.NonMinimalStrategy.Lomsac, values)


# ---------------------------------------------------------------
# MAGSAC
# ---------------------------------------------------------------

@unittest.skipUnless(homlib._HAS_MAGSAC, "MAGSAC not built")
class MAGSACSamplerTestCase(unittest.TestCase):
    def setUp(self):
        scene = _make_scene()
        self.x, self.y, self.scores, self.H, self.n_out = scene
        self.N = self.x.shape[1]
        self.solver = homlib.Wadenback3DV2026OneSided()
        self.sigma = 0.07

    # --- Samplers with scores ---

    @retry()
    def test_uniform_no_scores(self):
        mg = homlib.MAGSAC(
            self.solver, sigma_max=self.sigma)
        est, stats = mg.estimate(self.x, self.y)
        _assert_good(
            self, est, stats, self.H, self.N, self.n_out,
            check_inliers=False)

    @retry()
    def test_prosac_real_scores(self):
        mg = homlib.MAGSAC(
            self.solver, sigma_max=self.sigma,
            sampler=homlib.GCSamplerType.PROSAC)
        est, stats = mg.estimate(
            self.x, self.y, scores=self.scores)
        _assert_good(
            self, est, stats, self.H, self.N, self.n_out,
            check_inliers=False)

    @retry()
    def test_prosac_no_scores(self):
        mg = homlib.MAGSAC(
            self.solver, sigma_max=self.sigma,
            sampler=homlib.GCSamplerType.PROSAC)
        est, stats = mg.estimate(self.x, self.y)
        _assert_good(
            self, est, stats, self.H, self.N, self.n_out,
            check_inliers=False)

    @retry()
    def test_arsampler_real_scores(self):
        mg = homlib.MAGSAC(
            self.solver, sigma_max=self.sigma,
            sampler=homlib.GCSamplerType.ARSampler)
        est, stats = mg.estimate(
            self.x, self.y, scores=self.scores)
        _assert_good(
            self, est, stats, self.H, self.N, self.n_out,
            check_inliers=False)

    @retry()
    def test_importance_sampler_real_scores(self):
        mg = homlib.MAGSAC(
            self.solver, sigma_max=self.sigma,
            sampler=homlib.GCSamplerType.ImportanceSampler)
        est, stats = mg.estimate(
            self.x, self.y, scores=self.scores)
        _assert_good(
            self, est, stats, self.H, self.N, self.n_out,
            check_inliers=False)

    # --- Error on missing scores ---

    def test_arsampler_no_scores_raises(self):
        mg = homlib.MAGSAC(
            self.solver, sigma_max=self.sigma,
            sampler=homlib.GCSamplerType.ARSampler)
        with self.assertRaises(RuntimeError):
            mg.estimate(self.x, self.y)

    def test_importance_sampler_no_scores_raises(self):
        mg = homlib.MAGSAC(
            self.solver, sigma_max=self.sigma,
            sampler=homlib.GCSamplerType.ImportanceSampler)
        with self.assertRaises(RuntimeError):
            mg.estimate(self.x, self.y)

    # --- NonMinimalStrategy ---

    @retry()
    def test_non_minimal_lomsac(self):
        nms = homlib.NonMinimalStrategy.Lomsac
        mg = homlib.MAGSAC(
            self.solver, sigma_max=self.sigma,
            non_minimal_strategy=nms)
        est, stats = mg.estimate(self.x, self.y)
        _assert_good(
            self, est, stats, self.H, self.N, self.n_out,
            check_inliers=False)

    @retry()
    def test_non_minimal_multisample(self):
        nms = homlib.NonMinimalStrategy.MultiSample
        mg = homlib.MAGSAC(
            self.solver, sigma_max=self.sigma,
            non_minimal_strategy=nms)
        est, stats = mg.estimate(self.x, self.y)
        _assert_good(
            self, est, stats, self.H, self.N, self.n_out,
            check_inliers=False)

    # --- Settings ---

    def test_settings_defaults(self):
        s = homlib.MAGSACSettings()
        self.assertEqual(s.max_iterations, 5000)
        self.assertEqual(s.min_iterations, 50)
        self.assertAlmostEqual(s.confidence, 0.99)
        self.assertTrue(s.use_magsac_plus_plus)
        self.assertEqual(s.partition_number, 5)

    def test_settings_sampler_field(self):
        s = homlib.MAGSACSettings()
        s.sampler = homlib.GCSamplerType.ARSampler
        self.assertEqual(
            s.sampler, homlib.GCSamplerType.ARSampler)

    def test_settings_non_minimal_strategy_field(self):
        s = homlib.MAGSACSettings()
        nms = homlib.NonMinimalStrategy.MultiSample
        s.non_minimal_strategy = nms
        self.assertEqual(s.non_minimal_strategy, nms)

    # --- Namespace aliases ---

    def test_namespace_aliases(self):
        self.assertIs(
            homlib.MAGSAC.SamplerType,
            homlib.GCSamplerType)
        self.assertIs(
            homlib.MAGSAC.Settings,
            homlib.MAGSACSettings)
        self.assertIs(
            homlib.MAGSAC.Statistics,
            homlib.MAGSACStatistics)
        self.assertIs(
            homlib.MAGSAC.NonMinimalStrategy,
            homlib.NonMinimalStrategy)


# ---------------------------------------------------------------
# SupeRANSAC
# ---------------------------------------------------------------

if homlib._HAS_SUPERANSAC:
    SA_Sampler = homlib.SuperansacSamplerType
    SA_Score = homlib.SuperansacScoringType
    SA_LO = homlib.SuperansacLocalOptType


@unittest.skipUnless(homlib._HAS_SUPERANSAC, "SupeRANSAC not built")
class SuperansacSamplerTestCase(unittest.TestCase):
    def setUp(self):
        scene = _make_scene()
        self.x, self.y, self.scores, self.H, self.n_out = scene
        self.N = self.x.shape[1]
        self.solver = homlib.Wadenback3DV2026OneSided()
        self.thr = 0.07
        self.no_lo = dict(
            local_optimization=SA_LO.None_,
            final_optimization=SA_LO.None_,
        )

    # --- Samplers with scores ---

    @retry()
    def test_uniform_no_scores(self):
        sa = homlib.SupeRANSAC(
            self.solver, inlier_threshold=self.thr,
            sampler=SA_Sampler.Uniform,
            **self.no_lo)
        est, stats = sa.estimate(self.x, self.y)
        _assert_good(
            self, est, stats, self.H, self.N, self.n_out,
            check_inliers=False)

    @retry()
    def test_prosac_no_scores(self):
        sa = homlib.SupeRANSAC(
            self.solver, inlier_threshold=self.thr,
            sampler=SA_Sampler.PROSAC,
            **self.no_lo)
        est, stats = sa.estimate(self.x, self.y)
        _assert_good(
            self, est, stats, self.H, self.N, self.n_out,
            check_inliers=False)

    @retry()
    def test_prosac_with_scores(self):
        sa = homlib.SupeRANSAC(
            self.solver, inlier_threshold=self.thr,
            sampler=SA_Sampler.PROSAC,
            **self.no_lo)
        est, stats = sa.estimate(
            self.x, self.y, scores=self.scores)
        _assert_good(
            self, est, stats, self.H, self.N, self.n_out,
            check_inliers=False)

    def test_arsampler_with_scores(self):
        sa = homlib.SupeRANSAC(
            self.solver, inlier_threshold=self.thr,
            sampler=SA_Sampler.ARSampler,
            **self.no_lo)
        est, stats = sa.estimate(
            self.x, self.y, scores=self.scores)
        _assert_good(
            self, est, stats, self.H, self.N, self.n_out)

    def test_importance_sampler_with_scores(self):
        sa = homlib.SupeRANSAC(
            self.solver, inlier_threshold=self.thr,
            sampler=SA_Sampler.ImportanceSampler,
            **self.no_lo)
        est, stats = sa.estimate(
            self.x, self.y, scores=self.scores)
        _assert_good(
            self, est, stats, self.H, self.N, self.n_out,
            check_inliers=False)

    # --- Error on missing scores ---

    def test_arsampler_no_scores_raises(self):
        sa = homlib.SupeRANSAC(
            self.solver, inlier_threshold=self.thr,
            sampler=SA_Sampler.ARSampler,
            **self.no_lo)
        with self.assertRaises(RuntimeError):
            sa.estimate(self.x, self.y)

    def test_importance_sampler_no_scores_falls_back(self):
        sa = homlib.SupeRANSAC(
            self.solver, inlier_threshold=self.thr,
            sampler=SA_Sampler.ImportanceSampler,
            **self.no_lo)
        est, stats = sa.estimate(self.x, self.y)
        self.assertEqual(est.homography.shape, (3, 3))

    # --- Scoring types ---

    @retry()
    def test_scoring_msac(self):
        sa = homlib.SupeRANSAC(
            self.solver, inlier_threshold=self.thr,
            scoring=SA_Score.MSAC,
            **self.no_lo)
        est, stats = sa.estimate(self.x, self.y)
        _assert_good(
            self, est, stats, self.H, self.N, self.n_out,
            check_inliers=False)

    @retry()
    def test_scoring_ransac(self):
        sa = homlib.SupeRANSAC(
            self.solver, inlier_threshold=self.thr,
            scoring=SA_Score.RANSAC,
            **self.no_lo)
        est, stats = sa.estimate(self.x, self.y)
        _assert_good(
            self, est, stats, self.H, self.N, self.n_out,
            check_inliers=False)

    # --- Local optimization types ---

    @retry()
    def test_lo_irls(self):
        sa = homlib.SupeRANSAC(
            self.solver, inlier_threshold=self.thr,
            local_optimization=SA_LO.IRLS,
            final_optimization=SA_LO.None_)
        est, stats = sa.estimate(self.x, self.y)
        _assert_good(
            self, est, stats, self.H, self.N, self.n_out,
            check_inliers=False)

    @retry()
    def test_lo_lsq(self):
        sa = homlib.SupeRANSAC(
            self.solver, inlier_threshold=self.thr,
            local_optimization=SA_LO.LSQ,
            final_optimization=SA_LO.None_)
        est, stats = sa.estimate(self.x, self.y)
        _assert_good(
            self, est, stats, self.H, self.N, self.n_out,
            check_inliers=False)

    # --- Namespace aliases ---

    def test_namespace_aliases(self):
        SA = homlib.SupeRANSAC
        self.assertIs(SA.SamplerType, SA_Sampler)
        self.assertIs(
            SA.Settings, homlib.SuperansacSettings)
        self.assertIs(
            SA.Statistics, homlib.SuperansacStatistics)
        self.assertIs(SA.ScoringType, SA_Score)
        self.assertIs(SA.LocalOptType, SA_LO)

    def test_namespace_aliases_sub_settings(self):
        SA = homlib.SupeRANSAC
        self.assertIs(
            SA.ARSamplerSettings,
            homlib.SuperansacARSamplerSettings)
        self.assertIs(
            SA.LOSettings,
            homlib.SuperansacLOSettings)
        self.assertIs(
            SA.DegensacSettings,
            homlib.SuperansacDegensacSettings)

    # --- Enum completeness ---

    def test_sampler_enum_values(self):
        members = {e.name for e in SA_Sampler}
        expected = {
            "Uniform", "PROSAC", "NAPSAC",
            "ProgressiveNAPSAC",
            "ImportanceSampler", "ARSampler",
        }
        self.assertEqual(members, expected)

    def test_scoring_enum_values(self):
        members = {e.name for e in SA_Score}
        expected = {
            "RANSAC", "MSAC", "MAGSAC", "MINPRAN",
            "ACRANSAC", "GAU", "ML", "GRID",
        }
        self.assertEqual(members, expected)

    def test_local_opt_enum_values(self):
        members = {e.name for e in SA_LO}
        expected = {
            "None_", "LSQ", "IRLS", "NestedRANSAC",
            "GCRANSAC", "IteratedLMEDS", "CrossValidation",
        }
        self.assertEqual(members, expected)

    def test_termination_enum_values(self):
        members = {e.name for e in homlib.SuperansacTerminationType}
        self.assertIn("RANSAC", members)

    def test_inlier_selector_enum_values(self):
        members = {e.name for e in homlib.SuperansacInlierSelectorType}
        expected = {"None_", "SpacePartitioningRANSAC"}
        self.assertEqual(members, expected)

    def test_neighborhood_type_enum_values(self):
        members = {e.name for e in homlib.SuperansacNeighborhoodType}
        expected = {"Grid", "FLANN_KNN", "FLANN_Radius"}
        self.assertEqual(members, expected)

    # --- Settings: main fields ---

    def test_settings_defaults(self):
        s = homlib.SuperansacSettings()
        self.assertGreater(s.max_iterations, 0)
        self.assertGreater(s.min_iterations, 0)
        self.assertGreater(s.confidence, 0)
        self.assertGreater(s.inlier_threshold, 0)

    def test_settings_scalar_fields_writable(self):
        s = homlib.SuperansacSettings()
        s.inlier_threshold = 2.5
        self.assertAlmostEqual(s.inlier_threshold, 2.5)
        s.confidence = 0.95
        self.assertAlmostEqual(s.confidence, 0.95)
        s.min_iterations = 100
        self.assertEqual(s.min_iterations, 100)
        s.max_iterations = 10000
        self.assertEqual(s.max_iterations, 10000)
        s.top_k = 10
        self.assertEqual(s.top_k, 10)
        s.use_sprt = False
        self.assertFalse(s.use_sprt)

    def test_settings_bool_flags_writable(self):
        s = homlib.SuperansacSettings()
        s.local_optimization_inside_loop = False
        self.assertFalse(
            s.local_optimization_inside_loop)
        s.fundamental_orientation_check = True
        self.assertTrue(
            s.fundamental_orientation_check)
        s.homography_bundle_refinement = True
        self.assertTrue(
            s.homography_bundle_refinement)
        s.final_refinement_weighted = True
        self.assertTrue(s.final_refinement_weighted)

    def test_settings_refinement_fields_writable(self):
        s = homlib.SuperansacSettings()
        s.final_refinement_rounds = 5
        self.assertEqual(s.final_refinement_rounds, 5)
        s.final_refinement_schedule_start = 0.1
        self.assertAlmostEqual(
            s.final_refinement_schedule_start, 0.1)
        s.final_refinement_schedule_end = 0.9
        self.assertAlmostEqual(
            s.final_refinement_schedule_end, 0.9)
        s.homography_bundle_max_iterations = 50
        self.assertEqual(
            s.homography_bundle_max_iterations, 50)
        s.essential_cheirality_point_number = 3
        self.assertEqual(
            s.essential_cheirality_point_number, 3)

    def test_settings_enum_fields_writable(self):
        s = homlib.SuperansacSettings()
        s.sampler = SA_Sampler.PROSAC
        self.assertEqual(s.sampler, SA_Sampler.PROSAC)
        s.scoring = SA_Score.MSAC
        self.assertEqual(s.scoring, SA_Score.MSAC)
        s.local_optimization = SA_LO.IRLS
        self.assertEqual(
            s.local_optimization, SA_LO.IRLS)
        s.final_optimization = SA_LO.LSQ
        self.assertEqual(
            s.final_optimization, SA_LO.LSQ)

    # --- Settings: sub-settings structs ---

    def test_ar_sampler_settings_writable(self):
        s = homlib.SuperansacSettings()
        s.ar_sampler_settings.estimator_variance = 0.5
        self.assertAlmostEqual(
            s.ar_sampler_settings.estimator_variance,
            0.5)
        s.ar_sampler_settings.randomness = 0.02
        self.assertAlmostEqual(
            s.ar_sampler_settings.randomness, 0.02)

    def test_lo_settings_writable(self):
        s = homlib.SuperansacSettings()
        s.lo_settings.max_iterations = 100
        self.assertEqual(
            s.lo_settings.max_iterations, 100)
        s.lo_settings.graph_cut_number = 3
        self.assertEqual(
            s.lo_settings.graph_cut_number, 3)
        s.lo_settings.sample_size_multiplier = 7
        self.assertEqual(
            s.lo_settings.sample_size_multiplier, 7)
        s.lo_settings.spatial_coherence_weight = 0.14
        self.assertAlmostEqual(
            s.lo_settings.spatial_coherence_weight,
            0.14)

    def test_final_opt_settings_writable(self):
        s = homlib.SuperansacSettings()
        s.final_opt_settings.max_iterations = 200
        self.assertEqual(
            s.final_opt_settings.max_iterations, 200)

    def test_degensac_settings_writable(self):
        s = homlib.SuperansacSettings()
        s.degensac_settings.enabled = True
        self.assertTrue(s.degensac_settings.enabled)
        s.degensac_settings.homography_threshold = 2.0
        self.assertAlmostEqual(
            s.degensac_settings.homography_threshold,
            2.0)
        s.degensac_settings.inner_ransac_iterations = 50
        self.assertEqual(
            s.degensac_settings.inner_ransac_iterations,
            50)

    def test_neighborhood_settings_writable(self):
        s = homlib.SuperansacSettings()
        ns = s.neighborhood_settings
        ns.neighborhood_size = 10.0
        self.assertAlmostEqual(
            ns.neighborhood_size, 10.0)
        ns.neighborhood_grid_density = 8.0
        self.assertAlmostEqual(
            ns.neighborhood_grid_density, 8.0)
        ns.nearest_neighbor_number = 5
        self.assertEqual(
            ns.nearest_neighbor_number, 5)

    # --- Sub-settings constructors ---

    def test_ar_sampler_settings_standalone(self):
        ars = homlib.SuperansacARSamplerSettings()
        ars.estimator_variance = 1.0
        self.assertAlmostEqual(
            ars.estimator_variance, 1.0)

    def test_lo_settings_standalone(self):
        lo = homlib.SuperansacLOSettings()
        lo.max_iterations = 42
        self.assertEqual(lo.max_iterations, 42)

    def test_degensac_settings_standalone(self):
        d = homlib.SuperansacDegensacSettings()
        d.enabled = False
        self.assertFalse(d.enabled)

    def test_neighborhood_settings_standalone(self):
        n = homlib.SuperansacNeighborhoodSettings()
        n.neighborhood_grid_density = 6.0
        self.assertAlmostEqual(
            n.neighborhood_grid_density, 6.0)
