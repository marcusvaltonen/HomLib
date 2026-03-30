import unittest
from approvaltests import verify
from tests.helpers import verify_numpy_array

import numpy as np
import homlib


class FitzgibbonCVPR2001TestCase(unittest.TestCase):

    def test_estimate_fitzgibbon_cvpr_2001_two_sided_equal_results(self):
        """Check the expected result."""

        p1 = np.array([
            [-0.291910956401174,   0.998181103444196,   0.823771713189742,  -0.352037947619727,  -0.645473006258758],
            [0.827397315523355,   0.380888182504428,   0.834719771881603,  -0.327091812180983,  -0.378883379128205]
        ])

        p2 = np.array([
            [-0.380278821800649,   0.528335392478284,   0.957534928343294,   0.376937156586257,   0.462529355018110],
            [-0.340077827788294,  -0.680480654905349,   0.410935343442919,   0.882649958739618,   0.738720655679304]
        ])

        sols = homlib.estimate_fitzgibbon_cvpr_2001_two_sided_equal(p1, p2)
        self.assertEqual(len(sols), 1)
        self.assertAlmostEqual(sols[0].distortion_parameter, -0.6554778901775142, 12)
        self.assertAlmostEqual(sols[0].distortion_parameter2, 0.0, 12)
        # verify(verify_numpy_array(sols[0].homography))  # FIXME: Need to debug this error  RuntimeError: This machine has no reporter configuration

    def test_estimate_fitzgibbon_cvpr_2001_two_sided_equal_dimensions01(self):
        """Check that an exception is raised when dimensions are incorrect."""
        p1 = np.random.randn(2, 5)
        p2 = np.random.randn(3, 5)

        with self.assertRaises(TypeError):
            _ = homlib.estimate_fitzgibbon_cvpr_2001_two_sided_equal(p1, p2)

    def test_estimate_fitzgibbon_cvpr_2001_two_sided_equal_dimensions02(self):
        """Check that an exception is raised when dimensions are incorrect."""
        p1 = np.random.randn(2, 5)
        p2 = np.random.randn(2, 6)
        with self.assertRaises(ValueError):
            sols = homlib.estimate_fitzgibbon_cvpr_2001_two_sided_equal(p1, p2) # noqa
