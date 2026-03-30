import unittest
from approvaltests import verify
import pytest

import numpy as np
import homlib

from tests.helpers import verify_numpy_array


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

    @pytest.mark.skip(reason="Something strange with the Kukelova solver")
    def test_estimate_kukelova_cvpr_2015_two_sided_length(self):
        sols = homlib.estimate_kukelova_cvpr_2015_two_sided(self.p1, self.p2, False)
        self.assertEqual(len(sols), 5)

    @pytest.mark.skip(reason="Something strange with the Kukelova solver")
    def test_estimate_kukelova_cvpr_2015_two_sided_sol0(self):
        sols = homlib.estimate_kukelova_cvpr_2015_two_sided(self.p1, self.p2, False)
        self.assertEqual(len(sols), 5)
        self.assertAlmostEqual(sols[0].distortion_parameter, 0.2537509800067458, self.tol)
        self.assertAlmostEqual(sols[0].distortion_parameter2, -1.9785613160596929, self.tol)
        verify(verify_numpy_array(sols[0].homography))

    @pytest.mark.skip(reason="Something strange with the Kukelova solver")
    def test_estimate_kukelova_cvpr_2015_two_sided_sol1(self):
        sols = homlib.estimate_kukelova_cvpr_2015_two_sided(self.p1, self.p2, False)
        self.assertEqual(len(sols), 5)
        self.assertAlmostEqual(sols[1].distortion_parameter, -1.1832137508476386, self.tol)
        self.assertAlmostEqual(sols[1].distortion_parameter2, -1.8809663034629707, self.tol)
        verify(verify_numpy_array(sols[1].homography))

    @pytest.mark.skip(reason="Something strange with the Kukelova solver")
    def test_estimate_kukelova_cvpr_2015_two_sided_sol2(self):
        sols = homlib.estimate_kukelova_cvpr_2015_two_sided(self.p1, self.p2, False)
        self.assertEqual(len(sols), 5)
        self.assertAlmostEqual(sols[2].distortion_parameter, -1.7814746086320201, self.tol)
        self.assertAlmostEqual(sols[2].distortion_parameter2, -2.079301697963529, self.tol)
        verify(verify_numpy_array(sols[2].homography))

    @pytest.mark.skip(reason="Something strange with the Kukelova solver")
    def test_estimate_kukelova_cvpr_2015_two_sided_sol3(self):
        sols = homlib.estimate_kukelova_cvpr_2015_two_sided(self.p1, self.p2, False)
        self.assertEqual(len(sols), 5)
        self.assertAlmostEqual(sols[3].distortion_parameter2, -2.136402668559706, self.tol)
        self.assertAlmostEqual(sols[3].distortion_parameter2, 0.6928831549898077, self.tol)
        verify(verify_numpy_array(sols[3].homography))

    @pytest.mark.skip(reason="Something strange with the Kukelova solver")
    def test_estimate_kukelova_cvpr_2015_two_sided_sol4(self):
        sols = homlib.estimate_kukelova_cvpr_2015_two_sided(self.p1, self.p2, False)
        self.assertEqual(len(sols), 5)
        self.assertAlmostEqual(sols[4].distortion_parameter, -0.6554778901775545, self.tol)
        self.assertAlmostEqual(sols[4].distortion_parameter2, -0.6554778901775485, self.tol)
        verify(verify_numpy_array(sols[4].homography))

    def test_estimate_kukelova_cvpr_2015_two_sided_dimensions01(self):
        """Check that an exception is raised when dimensions are incorrect."""
        p1 = np.random.randn(2, 5)
        p2 = np.random.randn(3, 5)
        with self.assertRaises(TypeError):
            _ = homlib.estimate_kukelova_cvpr_2015_two_sided(p1, p2, False)

    def test_estimate_kukelova_cvpr_2015_two_sided_dimensions02(self):
        """Check that an exception is raised when dimensions are incorrect."""
        p1 = np.random.randn(2, 5)
        p2 = np.random.randn(2, 6)
        with self.assertRaises(ValueError):
            _ = homlib.estimate_kukelova_cvpr_2015_two_sided(p1, p2, False)
