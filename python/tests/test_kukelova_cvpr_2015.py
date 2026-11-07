import unittest

import numpy as np
import homlib


class KukelovaCVPR2015TestCase(unittest.TestCase):
    def setUp(self):
        self.p1 = np.array([
            [-0.291910956401174,   0.998181103444196,   0.823771713189742,  -0.352037947619727,  -0.645473006258758],
            [0.827397315523355,   0.380888182504428,   0.834719771881603,  -0.327091812180983,  -0.378883379128205]
        ])

        self.p2 = np.array([
            [-0.380278821800649,   0.528335392478284,   0.957534928343294,   0.376937156586257,   0.462529355018110],
            [-0.340077827788294,  -0.680480654905349,   0.410935343442919,   0.882649958739618,   0.738720655679304]
        ])

        self.tol = 10
        self.solver = homlib.KukelovaCVPR2015TwoSided()

    def test_estimate_kukelova_cvpr_2015_two_sided(self):
        sols = self.solver.estimate(self.p1, self.p2)
        self.assertEqual(len(sols), 5)

        expected_dp = [
            0.2537509800067458, -1.1832137508476386, -1.7814746086320201,
            -2.136402668559706, -0.6554778901775545,
        ]
        expected_dp2 = [
            -1.9785613160596929, -1.8809663034629707, -2.079301697963529,
            0.6928831549898077, -0.6554778901775485,
        ]
        expected_H = [
            np.array([[-1.84993859, 4.64801961, 0.72256438],
                      [-2.97548333, 4.41549851, 0.14438585],
                      [2.57097121, -6.03637522, -0.79542495]]),
            np.array([[-0.07243639, 0.3803624, 0.27070322],
                      [0.07137307, 0.35916405, 0.51141176],
                      [0.04851674, -0.51620915, -0.47045627]]),
            np.array([[0.16258481, 0.42614499, 0.33353189],
                      [0.3407417, 0.54233918, 0.50396602],
                      [-0.20919435, -0.52457951, -0.41537822]]),
            np.array([[0.34484818, 0.3490138, 0.35962541],
                      [0.57428832, 0.55155001, 0.5085264],
                      [0.99896668, 1.40954432, 1.14642718]]),
            np.array([[-0.0807097, -0.26843331, 0.13845705],
                      [0.20350208, -0.41058168, 0.57105471],
                      [-0.24702284, 0.11652806, 0.23191969]]),
        ]

        for i in range(5):
            self.assertAlmostEqual(sols[i].distortion_parameter, expected_dp[i], self.tol)
            self.assertAlmostEqual(sols[i].distortion_parameter2, expected_dp2[i], self.tol)
            np.testing.assert_allclose(sols[i].homography, expected_H[i], atol=1e-7)

    def test_estimate_kukelova_cvpr_2015_two_sided_dimensions01(self):
        """Check that an exception is raised when dimensions are incorrect."""
        p1 = np.random.randn(2, 5)
        p2 = np.random.randn(3, 5)
        with self.assertRaises(TypeError):
            _ = self.solver.estimate(p1, p2)

    def test_estimate_kukelova_cvpr_2015_two_sided_dimensions02(self):
        """Check that an exception is raised when dimensions are incorrect."""
        p1 = np.random.randn(2, 5)
        p2 = np.random.randn(2, 6)
        with self.assertRaises(ValueError):
            _ = self.solver.estimate(p1, p2)
