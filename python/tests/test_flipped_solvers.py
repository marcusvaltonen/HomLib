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


def assert_inverse(H_inv, H):
    """H_inv is the inverse of H, to the accuracy of an inverse: cond(H) * epsilon, relatively.

    Some solvers find nearly singular homographies (cond(H) up to 1e10 for Nakano's spurious roots):
    an absolute tolerance of 1e-10 is too tight for their inverses.
    """
    H_inv_ref = np.linalg.inv(H)
    tol = 100 * np.linalg.cond(H) * np.finfo(float).eps * np.abs(H_inv_ref).max()
    np.testing.assert_allclose(H_inv, H_inv_ref, rtol=0, atol=max(tol, 1e-10))


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


def verify_flipped(test, right_cls, left_cls, x1, x2):
    """The flipped solver returns the native solutions, with H inverted and k1, k2 swapped."""
    native_results = left_cls().estimate(x1, x2)
    flipped_results = right_cls().estimate(x2, x1)
    test.assertEqual(len(native_results), len(flipped_results))
    for n, f in zip(native_results, flipped_results):
        test.assertAlmostEqual(f.distortion_parameter, n.distortion_parameter2, places=10)
        test.assertAlmostEqual(f.distortion_parameter2, n.distortion_parameter, places=10)
        assert_inverse(f.homography, n.homography)


class FlippedSolverConsistencyTestCase(unittest.TestCase):
    """Flipped (left/default) vs native (right) solvers on valid samples."""

    def setUp(self):
        self.inst = make_instance(homlib.DistortionCase.ONE_SIDED_LEFT, 5)
        self.x1 = to_2n(self.inst.x1)
        self.x2 = to_2n(self.inst.x2)

    def _verify_consistency(self, right_cls, left_cls):
        native_results = left_cls().estimate(self.x1, self.x2)
        self.assertGreater(len(native_results), 0)
        self.assertTrue(has_gt_solution(native_results, self.inst))
        verify_flipped(self, right_cls, left_cls, self.x1, self.x2)

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


class FlippedNakanoNearlySingularTestCase(unittest.TestCase):
    """Samples on which Nakano's solver returns a nearly singular homography, a spurious root.

    Found among random samples (test_nakano): about 1 in 400 failed the count before the flipped
    solvers kept these solutions, and 1 in 2500 failed the comparison of the inverses before they
    were computed by LU and compared to the accuracy of an inverse.
    """

    def test_kept(self):
        # A spurious root with |det(H / |H|)| = 8e-15: the flipped solver used to drop it
        x1 = np.array([
            [-0.31245228179784473, -0.5439192567278298, -0.46297387620928454, -0.5056052181413829,
             -0.4379825191563144],
            [-0.5187985881007751, -0.21455288772488967, -0.6074444410777727, -0.470713350073152,
             -0.46084611686227783]])
        x2 = np.array([
            [-0.5210038135947531, -0.5819101636973875, -0.3934858030676815, -0.4469104836341385,
             -0.4879858431167529],
            [0.25473349904927123, -0.10898825546861696, 0.19240573494751254, 0.08815038522764963,
             0.12821040441977236]])
        self.assertEqual(len(homlib.NakanoICPR2025OneSided().estimate(x1, x2)), 3)
        verify_flipped(self, homlib.NakanoICPR2025OneSidedRight, homlib.NakanoICPR2025OneSided,
                       x1, x2)

    def test_inverse(self):
        # A spurious root with cond(H) = 1.4e6, whose inverse by cofactors was off by 1e-7
        x1 = np.array([
            [0.4322200772238785, 0.4964378638014161, 0.5434922462567146, 0.652388402240797,
             0.4032067401436877],
            [0.14739308768095294, -0.03145103528003382, 0.1391266466681737, -0.46596636268125574,
             0.21132589114262712]])
        x2 = np.array([
            [-0.5648037853299605, -0.5163723067925966, -0.44646415532699707, -0.42873558077531726,
             -0.5921543481813508],
            [0.5068915313499195, 0.2772239786330081, 0.45649840461268587, -0.12675039369708777,
             0.6046362952606202]])
        verify_flipped(self, homlib.NakanoICPR2025OneSidedRight, homlib.NakanoICPR2025OneSided,
                       x1, x2)


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
            assert_inverse(l.homography, r.homography)


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
            assert_inverse(l.homography, r.homography)


if __name__ == "__main__":
    unittest.main()
