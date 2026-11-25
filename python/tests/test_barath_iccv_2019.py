import unittest

import numpy as np
import homlib


def _points(inst):
    return np.array(inst.x1).T, np.array(inst.x2).T


def _ori(inst):
    return np.column_stack(inst.ori)


def _scale(inst):
    return np.column_stack(inst.scale)


class SIFTSolverTestCase(unittest.TestCase):
    """BarathICCV2019SIFT: 2 scale-and-orientation correspondences."""

    def test_minimal(self):
        config = homlib.ProblemConfig(homlib.DistortionCase.NO_DISTORTION, 0.0, 2)
        config.seed = 42
        inst = homlib.generate_problem_instance(config)
        x, y = _points(inst)
        ori = _ori(inst)
        scale = _scale(inst)
        solver = homlib.BarathICCV2019SIFT()
        sols = solver.estimate(x, y, ori=ori, scale=scale)
        self.assertGreaterEqual(len(sols), 1)
        self.assertLessEqual(len(sols), 2)
        best = min(sols, key=lambda s: inst.hom_error(s.homography))
        self.assertLess(inst.hom_error(best.homography), 1e-6)

    def test_feature_type(self):
        solver = homlib.BarathICCV2019SIFT()
        self.assertEqual(solver.feature_type, homlib.FeatureType.ScaleAndOrientation)

    def test_no_overdetermined(self):
        self.assertFalse(homlib.BarathICCV2019SIFT().supports_overdetermined)

    def test_multiple_instances(self):
        for seed in range(10):
            config = homlib.ProblemConfig(homlib.DistortionCase.NO_DISTORTION, 0.0, 2)
            config.seed = seed
            inst = homlib.generate_problem_instance(config)
            x, y = _points(inst)
            solver = homlib.BarathICCV2019SIFT()
            sols = solver.estimate(x, y, ori=_ori(inst), scale=_scale(inst))
            self.assertGreaterEqual(len(sols), 1, f"seed={seed}: no solutions")
            best = min(sols, key=lambda s: inst.hom_error(s.homography))
            self.assertLess(inst.hom_error(best.homography), 1e-6,
                            f"seed={seed}: hom_error={inst.hom_error(best.homography):.2e}")


class OrientationSolverTestCase(unittest.TestCase):
    """BarathICCV2019Orientation: 3 orientation correspondences."""

    def test_minimal(self):
        config = homlib.ProblemConfig(homlib.DistortionCase.NO_DISTORTION, 0.0, 3)
        config.seed = 42
        inst = homlib.generate_problem_instance(config)
        x, y = _points(inst)
        ori = _ori(inst)
        solver = homlib.BarathICCV2019Orientation()
        sols = solver.estimate(x, y, ori=ori)
        self.assertGreaterEqual(len(sols), 1)
        best = min(sols, key=lambda s: inst.hom_error(s.homography))
        self.assertLess(inst.hom_error(best.homography), 1e-6)

    def test_overdetermined(self):
        config = homlib.ProblemConfig(homlib.DistortionCase.NO_DISTORTION, 0.0, 10)
        inst = homlib.generate_problem_instance(config)
        x, y = _points(inst)
        ori = _ori(inst)
        solver = homlib.BarathICCV2019Orientation()
        sols = solver.estimate(x, y, ori=ori)
        self.assertEqual(len(sols), 1)
        self.assertLess(inst.hom_error(sols[0].homography), 1e-6)

    def test_feature_type(self):
        solver = homlib.BarathICCV2019Orientation()
        self.assertEqual(solver.feature_type, homlib.FeatureType.Orientation)

    def test_supports_overdetermined(self):
        self.assertTrue(homlib.BarathICCV2019Orientation().supports_overdetermined)

    def test_multiple_instances(self):
        for seed in range(10):
            config = homlib.ProblemConfig(homlib.DistortionCase.NO_DISTORTION, 0.0, 3)
            config.seed = seed
            inst = homlib.generate_problem_instance(config)
            x, y = _points(inst)
            solver = homlib.BarathICCV2019Orientation()
            sols = solver.estimate(x, y, ori=_ori(inst))
            self.assertGreaterEqual(len(sols), 1, f"seed={seed}: no solutions")
            best = min(sols, key=lambda s: inst.hom_error(s.homography))
            self.assertLess(inst.hom_error(best.homography), 1e-6,
                            f"seed={seed}: hom_error={inst.hom_error(best.homography):.2e}")


if __name__ == "__main__":
    unittest.main()
