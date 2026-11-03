import unittest

import numpy as np
import homlib


class Wadenback3DV2026OneSidedTestCase(unittest.TestCase):
    def setUp(self):
        self.p1 = np.array([
            [0.242104535827652,  0.319776775526576, -0.122936754370192,  0.303141152380151,  0.446905125282507],
            [-0.296297354299645, -0.114580331777873, -0.532655969828257, -0.295669384794041, -0.528127859697035]
        ])

        self.p2 = np.array([
            [-0.562787932824692, -0.612813481767545, -0.577247052139350, -0.544031097770073, -0.424602042497788],
            [-0.0769625591786919, -0.0311462459576149, -0.252034297773712, -0.0440402280844913,  0.0211895908990797]
        ])

        self.solver = homlib.Wadenback3DV2026OneSided(extra_check=False)

    def test_one_sided_length(self):
        sols = self.solver.estimate(self.p1, self.p2)
        self.assertEqual(len(sols), 3)

    def test_one_sided_distortion(self):
        sols = self.solver.estimate(self.p1, self.p2)
        self.assertAlmostEqual(sols[0].distortion_parameter, 0.0, 12)
        self.assertAlmostEqual(sols[0].distortion_parameter2, 5.11091, 4)
        self.assertAlmostEqual(sols[1].distortion_parameter, 0.0, 12)
        self.assertAlmostEqual(sols[1].distortion_parameter2, 2.73767, 4)
        self.assertAlmostEqual(sols[2].distortion_parameter, 0.0, 12)
        self.assertAlmostEqual(sols[2].distortion_parameter2, -0.3, 12)

    def test_one_sided_homography(self):
        sols = self.solver.estimate(self.p1, self.p2)
        expected_2 = np.array([
            [0.000227753625098076, -7.00643123791218e-05, -0.000525931983848800],
            [0.000429490139445578,  3.06112333501461e-05, -0.000156454396552730],
            [0.000133561563704977, -0.000486333541936380,  0.000545809522795618]
        ])
        np.testing.assert_allclose(sols[2].homography, expected_2, atol=1e-7)

    def test_one_sided_dimensions(self):
        p1 = np.random.randn(2, 5)
        p2 = np.random.randn(3, 5)
        with self.assertRaises(TypeError):
            self.solver.estimate(p1, p2)

    def test_one_sided_dimensions_cols(self):
        p1 = np.random.randn(2, 5)
        p2 = np.random.randn(2, 6)
        with self.assertRaises(ValueError):
            self.solver.estimate(p1, p2)


class Wadenback3DV2026TwoSidedEqualTestCase(unittest.TestCase):
    def setUp(self):
        self.p1 = np.array([
            [-0.093972209335993, -0.0969912533237537, -0.123299663772998, -0.300610529524079,  0.143634675965547],
            [0.0431909090094019,  0.0171160938161422,  0.500859415499538,  0.477995136756032,  0.402159621931429]
        ])

        self.p2 = np.array([
            [-0.0150485985127741, -0.0324017385815406, 0.335015888737872, 0.448124752741167, 0.0822796392160207],
            [0.504312576320796,   0.525321116499528,  0.182133325842435, 0.313117448200805, 0.0746306773026531]
        ])

        self.solver = homlib.Wadenback3DV2026TwoSidedEqual(extra_check=True)

    def test_two_sided_equal_length(self):
        sols = self.solver.estimate(self.p1, self.p2)
        self.assertEqual(len(sols), 1)

    def test_two_sided_equal_distortion(self):
        sols = self.solver.estimate(self.p1, self.p2)
        self.assertAlmostEqual(sols[0].distortion_parameter, -0.3, 12)
        self.assertAlmostEqual(sols[0].distortion_parameter2, -0.3, 12)

    def test_two_sided_equal_homography(self):
        sols = self.solver.estimate(self.p1, self.p2)
        expected = np.array([
            [-5.37070541971519e-06,  5.58776046545164e-06, -8.62056987828971e-07],
            [-5.32137209241095e-06, -4.87624482275744e-06,  3.51720726336518e-06],
            [1.80217660206107e-06,  2.74275008576624e-06,  7.02453977468668e-06]
        ])
        np.testing.assert_allclose(sols[0].homography, expected, atol=1e-7)


class Wadenback3DV2026TwoSidedTestCase(unittest.TestCase):
    def setUp(self):
        self.p1 = np.array([
            [-0.093972209335993, -0.0969912533237537, -0.123299663772998, -0.300610529524079,  0.143634675965547],
            [0.0431909090094019,  0.0171160938161422,  0.500859415499538,  0.477995136756032,  0.402159621931429]
        ])

        self.p2 = np.array([
            [-0.0150485985127741, -0.0324017385815406, 0.335015888737872, 0.448124752741167, 0.0822796392160207],
            [0.504312576320796,   0.525321116499528,  0.182133325842435, 0.313117448200805, 0.0746306773026531]
        ])

        self.solver = homlib.Wadenback3DV2026TwoSided(extra_check=True)

    def test_two_sided_length(self):
        sols = self.solver.estimate(self.p1, self.p2)
        self.assertEqual(len(sols), 1)

    def test_two_sided_distortion(self):
        sols = self.solver.estimate(self.p1, self.p2)
        self.assertAlmostEqual(sols[0].distortion_parameter, -0.3, 12)
        self.assertAlmostEqual(sols[0].distortion_parameter2, -0.3, 12)

    def test_two_sided_homography(self):
        sols = self.solver.estimate(self.p1, self.p2)
        expected = np.array([
            [-5.37070541971519e-06,  5.58776046545164e-06, -8.62056987828971e-07],
            [-5.32137209241095e-06, -4.87624482275744e-06,  3.51720726336518e-06],
            [1.80217660206107e-06,  2.74275008576624e-06,  7.02453977468668e-06]
        ])
        np.testing.assert_allclose(sols[0].homography, expected, atol=1e-7)
