import time
import unittest

import numpy as np
import homlib

from tests.helpers import normalize_H


@unittest.skipUnless(homlib._HAS_MAGSAC, "MAGSAC not built")
class MAGSACIntegrationTestCase(unittest.TestCase):
    def setUp(self):
        np.random.seed(42)

        self.H = np.array([
            [1.2, 0.1, 0.05],
            [0.05, 0.9, -0.02],
            [0.001, -0.001, 1.0]
        ])

        N = 100
        x = np.random.randn(2, N) * 0.3
        y_hom = self.H @ np.vstack((x, np.ones((1, N))))
        y = y_hom[:2] / y_hom[2]

        noise_std = 1e-4
        x += np.random.randn(2, N) * noise_std
        y += np.random.randn(2, N) * noise_std

        n_outliers = 20
        x[:, :n_outliers] = np.random.randn(2, n_outliers) * 10
        y[:, :n_outliers] = np.random.randn(2, n_outliers) * 10

        self.x = x
        self.y = y
        self.N = N
        self.n_outliers = n_outliers

    def test_magsac_wadenback_one_sided(self):
        solver = homlib.Wadenback3DV2026OneSided()
        mg = homlib.MAGSAC(solver, sigma_max=0.5)

        estimate, stats = mg.estimate(self.x, self.y)

        H_est = normalize_H(estimate.homography)
        H_gt = normalize_H(self.H)
        self.assertLess(np.linalg.norm(H_est - H_gt), 0.1)
        self.assertGreater(len(stats.inliers), self.N - self.n_outliers - 10)

    def test_magsac_fitzgibbon_one_sided(self):
        solver = homlib.FitzgibbonCVPR2001OneSided()
        mg = homlib.MAGSAC(solver, sigma_max=0.5)

        estimate, stats = mg.estimate(self.x, self.y)

        H_est = normalize_H(estimate.homography)
        H_gt = normalize_H(self.H)
        self.assertLess(np.linalg.norm(H_est - H_gt), 0.1)

    def test_magsac_settings_defaults(self):
        settings = homlib.MAGSACSettings()
        self.assertGreater(settings.sigma_max, 0)
        self.assertGreater(settings.confidence, 0)
        self.assertGreater(settings.max_iterations, 0)
        self.assertTrue(settings.use_magsac_plus_plus)

    def test_magsac_statistics_fields(self):
        solver = homlib.Wadenback3DV2026OneSided()
        mg = homlib.MAGSAC(solver, sigma_max=0.5)

        _, stats = mg.estimate(self.x, self.y)

        self.assertIsInstance(stats.inliers, list)
        self.assertGreater(stats.num_iterations, 0)

    def test_magsac_vs_lomsac(self):
        solver = homlib.Wadenback3DV2026OneSided()
        mg = homlib.MAGSAC(solver, sigma_max=0.5)
        lomsac = homlib.LOMSAC(solver, squared_inlier_threshold=0.005)

        H_gt = normalize_H(self.H)

        n_runs = 5
        mg_times, lomsac_times = [], []
        mg_inliers, lomsac_inliers = [], []

        for _ in range(n_runs):
            t0 = time.perf_counter()
            est_mg, stats_mg = mg.estimate(self.x, self.y)
            mg_times.append(time.perf_counter() - t0)
            mg_inliers.append(len(stats_mg.inliers))

            t0 = time.perf_counter()
            est_lo, stats_lo = lomsac.estimate(self.x, self.y)
            lomsac_times.append(time.perf_counter() - t0)
            lomsac_inliers.append(len(stats_lo.inlier_indices))

        H_mg = normalize_H(est_mg.homography)
        H_lo = normalize_H(est_lo.homography)
        self.assertLess(np.linalg.norm(H_mg - H_gt), 0.1)
        self.assertLess(np.linalg.norm(H_lo - H_gt), 0.1)

        print(f"\n  {'Backend':<14s} {'Avg time (ms)':>14s} {'Avg inliers':>12s}")
        print(f"  {'-'*42}")
        print(f"  {'MAGSAC':<14s} {np.mean(mg_times)*1000:14.2f} {np.mean(mg_inliers):12.1f}")
        print(f"  {'LOMSAC':<14s} {np.mean(lomsac_times)*1000:14.2f} {np.mean(lomsac_inliers):12.1f}")

    def test_magsac_plus_plus_vs_original(self):
        solver = homlib.Wadenback3DV2026OneSided()
        mgpp = homlib.MAGSAC(solver, sigma_max=0.5, use_magsac_plus_plus=True)
        mg_orig = homlib.MAGSAC(solver, sigma_max=0.5, use_magsac_plus_plus=False)

        est_pp, stats_pp = mgpp.estimate(self.x, self.y)
        est_orig, stats_orig = mg_orig.estimate(self.x, self.y)

        H_gt = normalize_H(self.H)
        H_pp = normalize_H(est_pp.homography)
        H_orig = normalize_H(est_orig.homography)

        self.assertLess(np.linalg.norm(H_pp - H_gt), 0.1)
        self.assertLess(np.linalg.norm(H_orig - H_gt), 0.5)


@unittest.skipIf(homlib._HAS_MAGSAC, "MAGSAC IS built — skip negative test")
class MAGSACImportErrorTestCase(unittest.TestCase):
    def test_magsac_import_error_without_build(self):
        solver = homlib.Wadenback3DV2026OneSided()
        with self.assertRaises(ImportError):
            homlib.MAGSAC(solver)
