import unittest

import homlib


class GenerateProblemInstanceTestCase(unittest.TestCase):
    def test_basic_generation(self):
        config = homlib.ProblemConfig(homlib.DistortionCase.ONE_SIDED_LEFT, 0.0, 5)
        instance = homlib.generate_problem_instance(config)
        self.assertEqual(len(instance.x1), 5)
        self.assertEqual(len(instance.x2), 5)

    def test_exact_ground_truth(self):
        config = homlib.ProblemConfig(homlib.DistortionCase.ONE_SIDED_LEFT, 0.0, 5)
        instance = homlib.generate_problem_instance(config)
        self.assertAlmostEqual(instance.hom_error(instance.posedata.homography), 0.0, 10)
        gt = instance.posedata
        self.assertAlmostEqual(
            instance.dist_error(gt.distortion_parameter, gt.distortion_parameter2), 0.0, 10
        )

    def test_noisy_ground_truth(self):
        config = homlib.ProblemConfig(homlib.DistortionCase.ONE_SIDED_LEFT, 1e-3, 10)
        instance = homlib.generate_problem_instance(config)
        self.assertEqual(len(instance.x1), 10)

    def test_distortion_cases_exist(self):
        for case in [
            homlib.DistortionCase.NO_DISTORTION,
            homlib.DistortionCase.ONE_SIDED_LEFT,
            homlib.DistortionCase.ONE_SIDED_RIGHT,
            homlib.DistortionCase.TWO_SIDED_EQUAL,
            homlib.DistortionCase.TWO_SIDED,
        ]:
            config = homlib.ProblemConfig(case, 0.0, 5)
            instance = homlib.generate_problem_instance(config)
            self.assertIsNotNone(instance)

    def test_problem_config_repr(self):
        config = homlib.ProblemConfig(homlib.DistortionCase.ONE_SIDED_LEFT, 0.0, 5)
        self.assertIn("ONE_SIDED_LEFT", repr(config))

    def test_problem_instance_repr(self):
        config = homlib.ProblemConfig(homlib.DistortionCase.ONE_SIDED_LEFT, 0.0, 5)
        instance = homlib.generate_problem_instance(config)
        self.assertIn("ProblemInstance", repr(instance))
