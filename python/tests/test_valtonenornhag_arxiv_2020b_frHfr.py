import unittest
from approvaltests import verify
from tests.helpers import verify_numpy_array

import numpy as np
import homlib


class ValtonenOrnhagArxiv2020BfrHfrTestCase(unittest.TestCase):
    def setUp(self):
        self.tol = 14

    def test_valtonenornhag_arxiv_2020b_frHfr(self):
        p1 = np.array([
            [-1.146281331283839,   1.050109134098990,   1.065996624259908],
            [-0.879951300627967,   0.620743795713172,   0.541580087112080],
        ])

        p2 = np.array([
            [0.663628650450811,   1.333268512835822,   1.318951998842419],
            [1.241359691717976,   0.068745345721370,   0.016786262835316],
        ])

        R1 = np.array([
            [-0.320761154096478,   0.935110133718446,   0.150603186685294],
            [-0.808554515336552,  -0.353152142517100,   0.470662469254195],
            [0.493307082608358,   0.029199350209406,   0.869365009760445],
        ])

        R2 = np.array([
            [0.420952545706761,  -0.744893719609945,  -0.517621773835547],
            [0.652124812245433,  -0.148125936460580,   0.743499788972085],
            [-0.630501533318403,  -0.650532130976893,   0.423409687005158],
        ])

        sols = homlib.get_valtonenornhag_arxiv_2020b_frHfr(
            np.asfortranarray(p1),
            np.asfortranarray(p2),
            np.asfortranarray(R1),
            np.asfortranarray(R2)
        )
        np.testing.assert_almost_equal(sols['lam'], -0.45054159850239145, self.tol)
        np.testing.assert_almost_equal(sols['f'], 5.756798219531202, self.tol)
        verify(verify_numpy_array(sols['H']))

    def test_valtonenornhag_arxiv_2020b_frHfr_dimensions01(self):
        """Check that an exception is raised when dimensions are incorrect."""
        p1 = np.random.randn(2, 3)
        p2 = np.random.randn(3, 3)
        R1 = np.random.randn(3, 3)
        R2 = np.random.randn(3, 3)
        with self.assertRaises(ValueError):
            sols = homlib.get_valtonenornhag_arxiv_2020b_frHfr( # noqa
                np.asfortranarray(p1),
                np.asfortranarray(p2),
                np.asfortranarray(R1),
                np.asfortranarray(R2)
            )

    def test_valtonenornhag_arxiv_2020b_frHfr_dimensions02(self):
        """Check that an exception is raised when dimensions are incorrect."""
        p1 = np.random.randn(2, 4)
        p2 = np.random.randn(2, 3)
        R1 = np.random.randn(3, 3)
        R2 = np.random.randn(3, 3)
        with self.assertRaises(ValueError):
            sols = homlib.get_valtonenornhag_arxiv_2020b_frHfr( # noqa
                np.asfortranarray(p1),
                np.asfortranarray(p2),
                np.asfortranarray(R1),
                np.asfortranarray(R2)
            )
