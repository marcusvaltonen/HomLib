import unittest
from approvaltests import verify
from tests.helpers import verify_numpy_array

import numpy as np
import homlib


class ValtonenOrnhagArxiv2020BfHfTestCase(unittest.TestCase):
    def setUp(self):
        self.tol = 14
        self.p1 = np.array([
            [2.003107199098924, -15.634084933471335],
            [-0.017087350257598, -7.041596829586987]
        ])

        self.p2 = np.array([
            [0.395688457559412, -0.012777594199286],
            [2.097270018093999, 0.988175585551782]
        ])

        self.R1 = np.array([
            [0.854451801803156, 0.080889542675225, 0.513194895026376],
            [0.251645807638113, 0.799754574299643, -0.545038538440134],
            [-0.454517882919365, 0.594852505057001, 0.662996222714662]
        ])

        self.R2 = np.array([
            [0.243935353955667, -0.895887070857591, -0.371324520279403],
            [0.945623784801441, 0.134783533871648, 0.296022054271079],
            [-0.215153900053711, -0.423343542843488, 0.880050591741408]
        ])

        self.sols = homlib.get_valtonenornhag_arxiv_2020b_fHf(
            np.asfortranarray(self.p1),
            np.asfortranarray(self.p2),
            np.asfortranarray(self.R1),
            np.asfortranarray(self.R2)
        )

    def test_valtonenornhag_arxiv_2020b_fHf_length(self):
        assert len(self.sols) == 4

    def test_valtonenornhag_arxiv_2020b_fHf_sol0(self):
        np.testing.assert_almost_equal(self.sols[0]['f'], -0.16963695303093723, self.tol)
        verify(verify_numpy_array(self.sols[0]['H']))

    def test_valtonenornhag_arxiv_2020b_fHf_sol1(self):
        np.testing.assert_almost_equal(self.sols[1]['f'], 0.34165415155423584, self.tol)
        verify(verify_numpy_array(self.sols[1]['H']))

    def test_valtonenornhag_arxiv_2020b_fHf_sol2(self):
        np.testing.assert_almost_equal(self.sols[2]['f'], 0.6535921559444265, self.tol)
        verify(verify_numpy_array(self.sols[2]['H']))

    def test_valtonenornhag_arxiv_2020b_fHf_sol3(self):
        np.testing.assert_almost_equal(self.sols[3]['f'], 3.3956095501687518, self.tol)
        verify(verify_numpy_array(self.sols[3]['H']))

    def test_valtonenornhag_arxiv_2020b_fHf_dimensions01(self):
        """Check that an exception is raised when dimensions are incorrect."""
        p1 = np.random.randn(2, 2)
        p2 = np.random.randn(3, 2)
        R1 = np.random.randn(3, 3)
        R2 = np.random.randn(3, 3)
        with self.assertRaises(ValueError):
            sols = homlib.get_valtonenornhag_arxiv_2020b_fHf( # noqa
                np.asfortranarray(p1),
                np.asfortranarray(p2),
                np.asfortranarray(R1),
                np.asfortranarray(R2)
            )

    def test_valtonenornhag_arxiv_2020b_fHf_dimensions02(self):
        """Check that an exception is raised when dimensions are incorrect."""
        p1 = np.random.randn(2, 3)
        p2 = np.random.randn(2, 2)
        R1 = np.random.randn(3, 3)
        R2 = np.random.randn(3, 3)
        with self.assertRaises(ValueError):
            sols = homlib.get_valtonenornhag_arxiv_2020b_fHf( # noqa
                np.asfortranarray(p1),
                np.asfortranarray(p2),
                np.asfortranarray(R1),
                np.asfortranarray(R2)
            )
