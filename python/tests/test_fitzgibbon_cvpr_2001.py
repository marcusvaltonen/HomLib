import unittest
from approvaltests import verify
from tests.helpers import verify_numpy_array

import numpy as np
import homlib


class FitzgibbonCVPR2001TestCase(unittest.TestCase):
    def setUp(self):
        self.p1 = np.array([
            [-0.291910956401174,   0.998181103444196,   0.823771713189742,  -0.352037947619727,  -0.645473006258758],
            [0.827397315523355,   0.380888182504428,   0.834719771881603,  -0.327091812180983,  -0.378883379128205]
        ])

        self.p2 = np.array([
            [-0.380278821800649,   0.528335392478284,   0.957534928343294,   0.376937156586257,   0.462529355018110],
            [-0.340077827788294,  -0.680480654905349,   0.410935343442919,   0.882649958739618,   0.738720655679304]
        ])

        self.tol = 15

    def test_fitzgibbon_cvpr_2001_results(self):
        """Check the expected result."""
        sols = homlib.get_fitzgibbon_cvpr_2001(np.asfortranarray(self.p1), np.asfortranarray(self.p2))
        np.testing.assert_almost_equal(sols['lam'], -0.6554778901775142, self.tol)
        verify(verify_numpy_array(sols['H']))

    def test_fitzgibbon_cvpr_2001_dimensions01(self):
        """Check that an exception is raised when dimensions are incorrect."""
        self.p1 = np.vstack((self.p1, np.ones((1, 5))))
        with self.assertRaises(ValueError):
            sols = homlib.get_fitzgibbon_cvpr_2001(np.asfortranarray(self.p1), np.asfortranarray(self.p2)) # noqa

    def test_fitzgibbon_cvpr_2001_dimensions02(self):
        """Check that an exception is raised when dimensions are incorrect."""
        self.p2 = np.vstack((self.p2, np.ones((1, 5))))
        with self.assertRaises(ValueError):
            sols = homlib.get_fitzgibbon_cvpr_2001(np.asfortranarray(self.p1), np.asfortranarray(self.p2)) # noqa
