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

    def test_lomsac_wadenback_two_sided_range(self):
        solver = homlib.Wadenback3DV2026TwoSided()
        self.assertIsNone(solver.k_range)

        # Range in the units of x and y, including the true distortion (0)
        solver.set_k_range(-0.1, 0.1)
        self.assertEqual(solver.k_range, (-0.1, 0.1))
        estimate, stats = homlib.LOMSAC(solver, squared_inlier_threshold=0.005).estimate(self.x, self.y)
        H_est = estimate.homography / estimate.homography[2, 2]
        self.assertLess(np.linalg.norm(H_est - self.H / self.H[2, 2]), 0.1)
        self.assertGreater(len(stats.inlier_indices), self.N - self.n_outliers - 10)

        # A range that excludes every minimal solution: no hypothesis, so no inliers
        solver.set_k_range(1e6, 2e6)
        _, stats = homlib.LOMSAC(solver, squared_inlier_threshold=0.005).estimate(self.x, self.y)
        self.assertEqual(len(stats.inlier_indices), 0)

        solver.reset_k_range()
        self.assertIsNone(solver.k_range)
        with self.assertRaises(ValueError):
            solver.set_k_range(0.1, -0.1)
        with self.assertRaises(ValueError):
            solver.set_k_range(k_min=float('nan'))
        with self.assertRaises(ValueError):
            solver.set_k_range(k_max=float('nan'))

    def test_lomsac_kukelova_two_sided_range(self):
        solver = homlib.KukelovaCVPR2015TwoSided()
        solver.set_k_range(-0.1, 0.1)
        estimate, stats = homlib.LOMSAC(solver, squared_inlier_threshold=0.005).estimate(self.x, self.y)
        H_est = estimate.homography / estimate.homography[2, 2]
        self.assertLess(np.linalg.norm(H_est - self.H / self.H[2, 2]), 0.1)
        self.assertGreater(len(stats.inlier_indices), self.N - self.n_outliers - 10)

        # A range that excludes every minimal solution: no hypothesis, so no inliers
        solver.set_k_range(k_min=1e6)
        _, stats = homlib.LOMSAC(solver, squared_inlier_threshold=0.005).estimate(self.x, self.y)
        self.assertEqual(len(stats.inlier_indices), 0)

    def test_lomsac_valtonenornhag_affine_range(self):
        # Noise-free affine-covariant problem, distorted in the first image
        config = homlib.ProblemConfig(homlib.DistortionCase.ONE_SIDED_RIGHT, 0.0, 50)
        instance = homlib.generate_problem_instance(config)
        x = np.array(instance.x1).T
        y = np.array(instance.x2).T
        A = np.array([np.asarray(a).flatten() for a in instance.A]).T
        k = instance.posedata.distortion_parameter

        # Range in the units of x and y, including the true distortion
        solver = homlib.ValtonenOrnhagICPR2026AffineSingleSided()
        solver.set_k_range(2 * k, k / 2)
        estimate, stats = homlib.LOMSAC(solver, squared_inlier_threshold=1e-6).estimate(x, y, A=A)
        self.assertAlmostEqual(estimate.distortion_parameter, k, places=6)
        self.assertGreater(len(stats.inlier_indices), 45)

        # A range that excludes every minimal solution: no hypothesis, so no inliers
        solver.set_k_range(k_min=1e6)
        _, stats = homlib.LOMSAC(solver, squared_inlier_threshold=1e-6).estimate(x, y, A=A)
        self.assertEqual(len(stats.inlier_indices), 0)

        # The direct estimate() applies the range too (this solver doesn't check it itself)
        solver.reset_k_range()
        self.assertGreater(len(solver.estimate(x[:, :2], y[:, :2], A=A[:, :2])), 0)
        solver.set_k_range(k_min=1e6)
        self.assertEqual(len(solver.estimate(x[:, :2], y[:, :2], A=A[:, :2])), 0)
        with self.assertRaises(TypeError):
            solver.set_k_range("a")

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
