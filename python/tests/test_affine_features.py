import unittest

import numpy as np

import homlib


class AffineFeaturesTestCase(unittest.TestCase):
    """The affine frames given to the affine solvers: each column of A is a row-major 2 x 2 matrix."""

    def test_row_major_affine_frames(self):
        solver = homlib.BarathVISAPP2016Affine()
        for _ in range(20):
            config = homlib.ProblemConfig(homlib.DistortionCase.NO_DISTORTION, 0.0, 2)
            instance = homlib.generate_problem_instance(config)
            x = np.array(instance.x1).T
            y = np.array(instance.x2).T
            A = np.array([np.asarray(a).flatten() for a in instance.A]).T   # row-major
            self.assertFalse(np.allclose(A[1], A[2]))   # not symmetric: the order matters
            poses = solver.estimate(x, y, A=A)
            self.assertEqual(len(poses), 1)
            self.assertLess(instance.hom_error(poses[0].homography), 1e-6)
            # The transposed frames give another homography
            poses_t = solver.estimate(x, y, A=A[[0, 2, 1, 3]])
            self.assertTrue(len(poses_t) == 0 or instance.hom_error(poses_t[0].homography) > 1e-3)


if __name__ == '__main__':
    unittest.main()
