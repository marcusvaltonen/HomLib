import unittest
from approvaltests import verify
from tests.helpers import verify_numpy_array

import numpy as np
import homlib


class ValtonenOrnhagArxiv2020ATestCase(unittest.TestCase):
    def setUp(self):
        self.p1 = np.array([
            [-0.170344517767220, 0.108864416575805, 0.661903907097742],
            [1.542599674500269, 0.956408496791784, 2.441464813306650]
        ])

        self.p2 = np.array([
            [0.729243664675995, 0.914374886276414, 3.761356659359217],
            [-0.571732334703832, -1.329303225626401, -1.507683466738961]
        ])

        self.R1 = np.array([
            [0.775811502326779, -0.525223981869245, -0.349651657692168],
            [0.591756227378766, 0.413366286756765, 0.692064216912979],
            [-0.218954516317697, -0.743819925682513, 0.631498881979806]
        ])

        self.R2 = np.array([
            [0.904829024395140, 0.344671896790600, 0.249971438718325],
            [-0.216541885342809, 0.878022117084280, -0.426833426295341],
            [-0.366597938488913, 0.332081986072128, 0.869095797954443]
        ])
        self.tol = 12

    def test_valtonenornhag_arxiv_2020a_fHf_results(self):
        """Check the expected result."""
        sols = homlib.get_valtonenornhag_arxiv_2020a_fHf(
            np.asfortranarray(self.p1),
            np.asfortranarray(self.p2),
            np.asfortranarray(self.R1),
            np.asfortranarray(self.R2)
        )
        np.testing.assert_almost_equal(sols['f'], 1.194848331672758, self.tol)
        verify(verify_numpy_array(sols['H']))

    def test_valtonenornhag_arxiv_2020a_fHf_dimensions01(self):
        """Check that an exception is raised when dimensions are incorrect."""
        self.p1 = np.vstack((self.p1, np.ones((1, 3))))
        with self.assertRaises(ValueError):
            sols = homlib.get_valtonenornhag_arxiv_2020a_fHf( # noqa
                np.asfortranarray(self.p1),
                np.asfortranarray(self.p2),
                np.asfortranarray(self.R1),
                np.asfortranarray(self.R2)
            )

    def test_valtonenornhag_arxiv_2020a_fHf_dimensions02(self):
        """Check that an exception is raised when dimensions are incorrect."""
        self.p2 = np.vstack((self.p2, np.ones((1, 3))))
        with self.assertRaises(ValueError):
            sols = homlib.get_valtonenornhag_arxiv_2020a_fHf( # noqa
                np.asfortranarray(self.p1),
                np.asfortranarray(self.p2),
                np.asfortranarray(self.R1),
                np.asfortranarray(self.R2)
            )
