import unittest

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

        solver = homlib.FitzgibbonCVPR2001TwoSidedEqual()
        sols = solver.estimate(p1, p2)
        self.assertEqual(len(sols), 1)
        self.assertAlmostEqual(sols[0].distortion_parameter, -0.6554778901775142, 12)
        self.assertAlmostEqual(sols[0].distortion_parameter2, -0.6554778901775142, 12)
        expected = np.array([
            [0.07731617, 0.25714674, -0.13263548],
            [-0.19494561, 0.39331832, -0.54704409],
            [0.23663649, -0.11162851, -0.22216837]
        ])
        np.testing.assert_allclose(sols[0].homography, expected, atol=1e-7)

    def test_estimate_fitzgibbon_cvpr_2001_two_sided_equal_dimensions01(self):
        """Check that an exception is raised when dimensions are incorrect."""
        p1 = np.random.randn(2, 5)
        p2 = np.random.randn(3, 5)

        solver = homlib.FitzgibbonCVPR2001TwoSidedEqual()
        with self.assertRaises(TypeError):
            _ = solver.estimate(p1, p2)

    def test_estimate_fitzgibbon_cvpr_2001_two_sided_equal_dimensions02(self):
        """Check that an exception is raised when dimensions are incorrect."""
        p1 = np.random.randn(2, 5)
        p2 = np.random.randn(2, 6)
        solver = homlib.FitzgibbonCVPR2001TwoSidedEqual()
        with self.assertRaises(ValueError):
            sols = solver.estimate(p1, p2) # noqa


class FitzgibbonOverdeterminedTestCase(unittest.TestCase):
    """Fitzgibbon solvers with N > 5 (overdetermined)."""

    @staticmethod
    def _points(inst):
        return np.array(inst.x1).T, np.array(inst.x2).T

    def test_one_sided_overdetermined(self):
        config = homlib.ProblemConfig(homlib.DistortionCase.ONE_SIDED_LEFT, 0.0, 20)
        inst = homlib.generate_problem_instance(config)
        x, y = self._points(inst)
        solver = homlib.FitzgibbonCVPR2001OneSided()
        sols = solver.estimate(x, y)
        self.assertEqual(len(sols), 1)
        self.assertAlmostEqual(
            sols[0].distortion_parameter2, inst.posedata.distortion_parameter2, places=2)

    def test_one_sided_fewer_than_5(self):
        solver = homlib.FitzgibbonCVPR2001OneSided()
        x = np.random.randn(2, 4)
        sols = solver.estimate(x, x)
        self.assertEqual(len(sols), 0)

    def test_two_sided_equal_overdetermined(self):
        config = homlib.ProblemConfig(homlib.DistortionCase.TWO_SIDED_EQUAL, 0.0, 20)
        inst = homlib.generate_problem_instance(config)
        x, y = self._points(inst)
        solver = homlib.FitzgibbonCVPR2001TwoSidedEqual()
        sols = solver.estimate(x, y)
        self.assertEqual(len(sols), 1)
        self.assertAlmostEqual(
            sols[0].distortion_parameter, inst.posedata.distortion_parameter, places=2)

    def test_minimal_case_unchanged(self):
        """Regression: minimal N=5 still produces the same result as before."""
        config = homlib.ProblemConfig(homlib.DistortionCase.ONE_SIDED_LEFT, 0.0, 5)
        inst = homlib.generate_problem_instance(config)
        x, y = self._points(inst)
        solver = homlib.FitzgibbonCVPR2001OneSided()
        sols = solver.estimate(x, y)
        self.assertGreaterEqual(len(sols), 1)
        self.assertAlmostEqual(
            sols[0].distortion_parameter2, inst.posedata.distortion_parameter2, places=3)
