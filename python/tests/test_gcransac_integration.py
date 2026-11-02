import time
import unittest

import numpy as np
import homlib


@unittest.skipUnless(homlib._HAS_GCRANSAC, "GC-RANSAC not built")
class GCRANSACIntegrationTestCase(unittest.TestCase):
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

    def test_gcransac_wadenback_one_sided(self):
        solver = homlib.Wadenback3DV2026OneSided()
        gc = homlib.GCRANSAC(solver, threshold=0.07)

        estimate, stats = gc.estimate(self.x, self.y)

        H_est = estimate.homography / estimate.homography[2, 2]
        H_gt = self.H / self.H[2, 2]
        self.assertLess(np.linalg.norm(H_est - H_gt), 0.1)
        self.assertGreater(len(stats.inliers), self.N - self.n_outliers - 10)

    def test_gcransac_fitzgibbon_one_sided(self):
        solver = homlib.FitzgibbonCVPR2001OneSided()
        gc = homlib.GCRANSAC(solver, threshold=0.07)

        estimate, stats = gc.estimate(self.x, self.y)

        H_est = estimate.homography / estimate.homography[2, 2]
        H_gt = self.H / self.H[2, 2]
        self.assertLess(np.linalg.norm(H_est - H_gt), 0.1)

    def test_gcransac_settings_defaults(self):
        settings = homlib.GCRANSACSettings()
        self.assertGreater(settings.threshold, 0)
        self.assertGreater(settings.confidence, 0)
        self.assertGreater(settings.max_iterations, 0)

    def test_gcransac_statistics_fields(self):
        solver = homlib.Wadenback3DV2026OneSided()
        gc = homlib.GCRANSAC(solver, threshold=0.07)

        _, stats = gc.estimate(self.x, self.y)

        self.assertIsInstance(stats.inliers, list)
        self.assertGreater(stats.num_iterations, 0)

    def test_gcransac_vs_lomsac(self):
        solver = homlib.Wadenback3DV2026OneSided()
        gc = homlib.GCRANSAC(solver, threshold=0.07)
        lomsac = homlib.LOMSAC(solver, squared_inlier_threshold=0.005)

        H_gt = self.H / self.H[2, 2]

        n_runs = 5
        gc_times, lomsac_times = [], []
        gc_inliers, lomsac_inliers = [], []

        for _ in range(n_runs):
            t0 = time.perf_counter()
            est_gc, stats_gc = gc.estimate(self.x, self.y)
            gc_times.append(time.perf_counter() - t0)
            gc_inliers.append(len(stats_gc.inliers))

            t0 = time.perf_counter()
            est_lo, stats_lo = lomsac.estimate(self.x, self.y)
            lomsac_times.append(time.perf_counter() - t0)
            lomsac_inliers.append(len(stats_lo.inlier_indices))

        H_gc = est_gc.homography / est_gc.homography[2, 2]
        H_lo = est_lo.homography / est_lo.homography[2, 2]
        self.assertLess(np.linalg.norm(H_gc - H_gt), 0.1)
        self.assertLess(np.linalg.norm(H_lo - H_gt), 0.1)

        print(f"\n  {'Backend':<14s} {'Avg time (ms)':>14s} {'Avg inliers':>12s}")
        print(f"  {'-'*42}")
        print(f"  {'GC-RANSAC':<14s} {np.mean(gc_times)*1000:14.2f} {np.mean(gc_inliers):12.1f}")
        print(f"  {'LOMSAC':<14s} {np.mean(lomsac_times)*1000:14.2f} {np.mean(lomsac_inliers):12.1f}")


@unittest.skipIf(homlib._HAS_GCRANSAC, "GC-RANSAC IS built — skip negative test")
class GCRANSACImportErrorTestCase(unittest.TestCase):
    def test_gcransac_import_error_without_build(self):
        solver = homlib.Wadenback3DV2026OneSided()
        with self.assertRaises(ImportError):
            homlib.GCRANSAC(solver)
