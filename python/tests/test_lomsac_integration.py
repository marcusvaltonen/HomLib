import unittest

import numpy as np
import homlib


class LomsacIntegrationTestCase(unittest.TestCase):
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

    def test_lomsac_wadenback_one_sided(self):
        solver = homlib.Wadenback3DV2026OneSided()
        lomsac = homlib.LOMSAC(solver, squared_inlier_threshold=0.005)

        estimate, stats = lomsac.estimate(self.x, self.y)

        H_est = estimate.homography / estimate.homography[2, 2]
        H_gt = self.H / self.H[2, 2]
        self.assertLess(np.linalg.norm(H_est - H_gt), 0.1)

        self.assertGreater(len(stats.inlier_indices), self.N - self.n_outliers - 10)

    def test_lomsac_fitzgibbon_one_sided(self):
        solver = homlib.FitzgibbonCVPR2001OneSided()
        lomsac = homlib.LOMSAC(solver, squared_inlier_threshold=0.005)

        estimate, stats = lomsac.estimate(self.x, self.y)

        H_est = estimate.homography / estimate.homography[2, 2]
        H_gt = self.H / self.H[2, 2]
        self.assertLess(np.linalg.norm(H_est - H_gt), 0.1)

    def test_lomsac_options_defaults(self):
        options = homlib.LORansacOptions()
        self.assertGreater(options.max_num_iterations, 0)
        self.assertGreater(options.success_probability, 0)

    def test_ransac_statistics_fields(self):
        solver = homlib.Wadenback3DV2026OneSided()
        lomsac = homlib.LOMSAC(solver, squared_inlier_threshold=0.005)

        _, stats = lomsac.estimate(self.x, self.y)

        self.assertIsInstance(stats.inlier_indices, list)
        self.assertGreater(stats.num_iterations, 0)
