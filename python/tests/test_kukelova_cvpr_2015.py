import unittest
from approvaltests import verify
from tests.helpers import verify_numpy_array

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

    def test_get_kukelova_cvpr_2015_length(self):
        sols = homlib.get_kukelova_cvpr_2015(np.asfortranarray(self.p1), np.asfortranarray(self.p2))
        assert len(sols) == 5

    def test_get_kukelova_cvpr_2015_sol0(self):
        sols = homlib.get_kukelova_cvpr_2015(np.asfortranarray(self.p1), np.asfortranarray(self.p2))
        np.testing.assert_almost_equal(sols[0]['lam1'], 0.2537509800067458, self.tol)
        np.testing.assert_almost_equal(sols[0]['lam2'], -1.9785613160596929, self.tol)
        verify(verify_numpy_array(sols[0]['H']))

    def test_get_kukelova_cvpr_2015_sol1(self):
        sols = homlib.get_kukelova_cvpr_2015(np.asfortranarray(self.p1), np.asfortranarray(self.p2))
        np.testing.assert_almost_equal(sols[1]['lam1'], -1.1832137508476386, self.tol)
        np.testing.assert_almost_equal(sols[1]['lam2'], -1.8809663034629707, self.tol)
        verify(verify_numpy_array(sols[1]['H']))

    def test_get_kukelova_cvpr_2015_sol2(self):
        sols = homlib.get_kukelova_cvpr_2015(np.asfortranarray(self.p1), np.asfortranarray(self.p2))
        np.testing.assert_almost_equal(sols[2]['lam1'], -1.7814746086320201, self.tol)
        np.testing.assert_almost_equal(sols[2]['lam2'], -2.079301697963529, self.tol)
        verify(verify_numpy_array(sols[2]['H']))

    def test_get_kukelova_cvpr_2015_sol3(self):
        sols = homlib.get_kukelova_cvpr_2015(np.asfortranarray(self.p1), np.asfortranarray(self.p2))
        np.testing.assert_almost_equal(sols[3]['lam1'], -2.136402668559706, self.tol)
        np.testing.assert_almost_equal(sols[3]['lam2'], 0.6928831549898077, self.tol)
        verify(verify_numpy_array(sols[3]['H']))

    def test_get_kukelova_cvpr_2015_sol4(self):
        sols = homlib.get_kukelova_cvpr_2015(np.asfortranarray(self.p1), np.asfortranarray(self.p2))
        np.testing.assert_almost_equal(sols[4]['lam1'], -0.6554778901775545, self.tol)
        np.testing.assert_almost_equal(sols[4]['lam2'], -0.6554778901775485, self.tol)
        verify(verify_numpy_array(sols[4]['H']))

    def test_kukelova_cvpr_2015_dimensions01(self):
        """Check that an exception is raised when dimensions are incorrect."""
        self.p1 = np.vstack((self.p1, np.ones((1, 5))))
        with self.assertRaises(ValueError):
            sols = homlib.get_kukelova_cvpr_2015(np.asfortranarray(self.p1), np.asfortranarray(self.p2)) # noqa

    def test_kukelova_cvpr_2015_dimensions02(self):
        """Check that an exception is raised when dimensions are incorrect."""
        self.p2 = np.vstack((self.p2, np.ones((1, 5))))
        with self.assertRaises(ValueError):
            sols = homlib.get_kukelova_cvpr_2015(np.asfortranarray(self.p1), np.asfortranarray(self.p2)) # noqa
