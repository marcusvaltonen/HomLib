import unittest

import numpy as np
import homlib


def _points(inst):
    return np.array(inst.x1).T, np.array(inst.x2).T


def _affine(inst):
    return np.array([np.asarray(a).flatten() for a in inst.A]).T


def _ori(inst):
    return np.column_stack(inst.ori)


class AffineOverdeterminedTestCase(unittest.TestCase):

    def test_overdetermined(self):
        config = homlib.ProblemConfig(homlib.DistortionCase.ONE_SIDED_RIGHT, 0.0, 10)
        inst = homlib.generate_problem_instance(config)
        x, y = _points(inst)
        A = _affine(inst)
        solver = homlib.ValtonenOrnhagICPR2026AffineOneSidedRight()
        sols = solver.estimate(x, y, A=A)
        self.assertEqual(len(sols), 1)
        self.assertAlmostEqual(
            sols[0].distortion_parameter, inst.posedata.distortion_parameter, places=2)

    def test_minimal_regression(self):
        config = homlib.ProblemConfig(homlib.DistortionCase.ONE_SIDED_RIGHT, 0.0, 2)
        inst = homlib.generate_problem_instance(config)
        x, y = _points(inst)
        A = _affine(inst)
        solver = homlib.ValtonenOrnhagICPR2026AffineOneSidedRight()
        sols = solver.estimate(x, y, A=A)
        self.assertGreaterEqual(len(sols), 1)
        self.assertAlmostEqual(
            sols[0].distortion_parameter, inst.posedata.distortion_parameter, places=3)


class OrientationOverdeterminedTestCase(unittest.TestCase):

    def test_overdetermined(self):
        config = homlib.ProblemConfig(homlib.DistortionCase.ONE_SIDED_RIGHT, 0.0, 15)
        inst = homlib.generate_problem_instance(config)
        x, y = _points(inst)
        ori = _ori(inst)
        solver = homlib.ValtonenOrnhagICPR2026OrientationOneSidedRight()
        sols = solver.estimate(x, y, ori=ori)
        self.assertEqual(len(sols), 1)
        self.assertAlmostEqual(
            sols[0].distortion_parameter, inst.posedata.distortion_parameter, places=2)

    def test_minimal_regression(self):
        # Three correspondences, the minimal sample: up to three solutions, one of them the true one
        config = homlib.ProblemConfig(homlib.DistortionCase.ONE_SIDED_RIGHT, 0.0, 3)
        config.seed = 42
        inst = homlib.generate_problem_instance(config)
        x, y = _points(inst)
        ori = _ori(inst)
        solver = homlib.ValtonenOrnhagICPR2026OrientationOneSidedRight()
        sols = solver.estimate(x, y, ori=ori)
        self.assertGreaterEqual(len(sols), 1)
        best = min(sols, key=lambda s: abs(s.distortion_parameter - inst.posedata.distortion_parameter))
        self.assertLess(abs(best.distortion_parameter - inst.posedata.distortion_parameter), 1e-6)
        self.assertLess(inst.hom_error(best.homography), 1e-6)   # H maps x to y, not y to x
