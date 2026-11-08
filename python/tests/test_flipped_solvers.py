import unittest

import numpy as np
import homlib


def make_instance(distortion, n_points):
    config = homlib.ProblemConfig(distortion, 0.0, n_points)
    return homlib.generate_problem_instance(config)


def to_2n(pts):
    """Convert list of (2,) arrays to [2, n] numpy array."""
    return np.column_stack(pts)


def to_A_4n(A_list):
    """Convert list of (2,2) affine matrices to [4, n] numpy array."""
    return np.column_stack([a.flatten() for a in A_list])


def has_gt_solution(results, inst, hom_tol=1e-6, dist_tol=1e-6):
    gt = inst.posedata
    H_gt = gt.homography / np.linalg.norm(gt.homography)
    if H_gt[2, 2] < 0:
        H_gt *= -1.0
    for r in results:
        H_est = r.homography / np.linalg.norm(r.homography)
        if H_est[2, 2] < 0:
            H_est *= -1.0
        h_err = np.linalg.norm(H_gt - H_est)
        d_err = 0.5 * (
            abs(gt.distortion_parameter - r.distortion_parameter)
            + abs(gt.distortion_parameter2 - r.distortion_parameter2)
        )
        if h_err < hom_tol and d_err < dist_tol:
            return True
    return False


class FlippedSolverConsistencyTestCase(unittest.TestCase):
    """Flipped (left/default) vs native (right) solvers on valid samples."""

    def setUp(self):
        self.inst = make_instance(homlib.DistortionCase.ONE_SIDED_LEFT, 5)
        self.x1 = to_2n(self.inst.x1)
        self.x2 = to_2n(self.inst.x2)

    def _verify_consistency(self, right_cls, left_cls):
        native_results = left_cls().estimate(self.x1, self.x2)
        flipped_results = right_cls().estimate(self.x2, self.x1)

        self.assertGreater(len(native_results), 0)
        self.assertTrue(has_gt_solution(native_results, self.inst))

        self.assertEqual(len(native_results), len(flipped_results))

        for n, f in zip(native_results, flipped_results):
            self.assertAlmostEqual(
                f.distortion_parameter,
                n.distortion_parameter2, places=10,
            )
            self.assertAlmostEqual(
                f.distortion_parameter2,
                n.distortion_parameter, places=10,
            )
            H_inv = np.linalg.inv(n.homography)
            np.testing.assert_allclose(
                f.homography, H_inv, atol=1e-10,
            )

    def test_fitzgibbon(self):
        self._verify_consistency(
            homlib.FitzgibbonCVPR2001OneSidedRight,
            homlib.FitzgibbonCVPR2001OneSided,
        )

    def test_nakano(self):
        self._verify_consistency(
            homlib.NakanoICPR2025OneSidedRight,
            homlib.NakanoICPR2025OneSided,
        )

    def test_wadenback(self):
        self._verify_consistency(
            homlib.Wadenback3DV2026OneSidedRight,
            homlib.Wadenback3DV2026OneSided,
        )


class FlippedAffineConsistencyTestCase(unittest.TestCase):
    """Flipped affine solver: left vs native right on valid samples."""

    def test_affine(self):
        inst = make_instance(homlib.DistortionCase.ONE_SIDED_RIGHT, 2)
        x1 = to_2n(inst.x1)
        x2 = to_2n(inst.x2)
        A = to_A_4n(inst.A)

        A_inv = np.column_stack([
            np.linalg.inv(a).flatten() for a in inst.A
        ])

        right = homlib.ValtonenOrnhagICPR2026AffineOneSidedRight()
        rr = right.estimate(x1, x2, A=A)
        left = homlib.ValtonenOrnhagICPR2026AffineOneSided()
        lr = left.estimate(x2, x1, A=A_inv)

        self.assertGreater(len(rr), 0)
        self.assertTrue(has_gt_solution(rr, inst))

        self.assertEqual(len(rr), len(lr))
        for r, l in zip(rr, lr):
            self.assertAlmostEqual(
                l.distortion_parameter,
                r.distortion_parameter2, places=10,
            )
            self.assertAlmostEqual(
                l.distortion_parameter2,
                r.distortion_parameter, places=10,
            )
            np.testing.assert_allclose(
                l.homography,
                np.linalg.inv(r.homography), atol=1e-10,
            )


class FlippedOrientationConsistencyTestCase(unittest.TestCase):
    """Flipped orientation solver: left vs native right on valid samples."""

    def test_orientation(self):
        inst = make_instance(homlib.DistortionCase.ONE_SIDED_RIGHT, 4)
        x1 = to_2n(inst.x1)
        x2 = to_2n(inst.x2)
        ori = to_2n(inst.ori)

        ori_swapped = np.array([ori[1], ori[0]])

        right = homlib.ValtonenOrnhagICPR2026OrientationOneSidedRight()
        rr = right.estimate(x1, x2, ori=ori)
        left = homlib.ValtonenOrnhagICPR2026OrientationOneSided()
        lr = left.estimate(x2, x1, ori=ori_swapped)

        self.assertGreater(len(rr), 0)
        self.assertTrue(has_gt_solution(rr, inst))

        self.assertEqual(len(rr), len(lr))
        for r, l in zip(rr, lr):
            self.assertAlmostEqual(
                l.distortion_parameter,
                r.distortion_parameter2, places=10,
            )
            self.assertAlmostEqual(
                l.distortion_parameter2,
                r.distortion_parameter, places=10,
            )
            np.testing.assert_allclose(
                l.homography,
                np.linalg.inv(r.homography), atol=1e-10,
            )


if __name__ == "__main__":
    unittest.main()
