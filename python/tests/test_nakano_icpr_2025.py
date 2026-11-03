import unittest

import numpy as np
import homlib


class NakanoICPR2025TestCase(unittest.TestCase):
    def setUp(self):
        self.p1 = np.array([
            [-0.514683699631253, -0.216333891824487,  0.0149270190286551, -0.299909368807735, -0.0667164585547235],
            [-0.0622906242414419, 0.236104964783804, -0.466167744934018,  -0.631181204140851,  0.113775620793983]
        ])

        self.p2 = np.array([
            [-0.173489966290854, -0.276098980665827,  0.389704168705386,  0.307091630022996, -0.0999459278021077],
            [0.161874840375985,  0.528540478431576,  0.291358438524245,  0.0216079362297207, 0.583861165740761]
        ])

        self.solver = homlib.NakanoICPR2025OneSided(extra_check=False)

    def test_length(self):
        sols = self.solver.estimate(self.p1, self.p2)
        self.assertEqual(len(sols), 1)

    def test_distortion(self):
        sols = self.solver.estimate(self.p1, self.p2)
        self.assertAlmostEqual(sols[0].distortion_parameter, 0.0, 12)
        self.assertAlmostEqual(sols[0].distortion_parameter2, -0.3, 12)

    def test_homography(self):
        sols = self.solver.estimate(self.p1, self.p2)
        expected = np.array([
            [-0.966294723999024,  1.424868744982000, -0.0570986060219543],
            [-1.118034960806880, -0.916112342444385, -0.9604495174992110],
            [0.851682151867801,  0.527896230136481, -1.5204980422764300]
        ])
        np.testing.assert_allclose(sols[0].homography, expected, atol=1e-7)

    def test_dimensions(self):
        p1 = np.random.randn(2, 5)
        p2 = np.random.randn(3, 5)
        with self.assertRaises(TypeError):
            self.solver.estimate(p1, p2)

    def test_dimensions_cols(self):
        p1 = np.random.randn(2, 5)
        p2 = np.random.randn(2, 6)
        with self.assertRaises(ValueError):
            self.solver.estimate(p1, p2)
