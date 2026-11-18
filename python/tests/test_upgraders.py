"""Tests for feature upgraders: AffineDecompositionUpgrader and BarathAffineUpgrader.

Uses homlib.generate_problem_instance() as ground truth: it outputs matched
(A, ori) pairs where A is the relative affine and ori = [theta1, theta2] from
the affine2sift decomposition, with q = sqrt(det(A)) as the relative scale.
"""

import sys
import unittest
from pathlib import Path

import numpy as np

import homlib

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "benchmark"))
from common import AffineMatchResult, ScaleOrientationMatchResult  # noqa: E402
from features import AffineDecompositionUpgrader, BarathAffineUpgrader  # noqa: E402


def _rotation_matrix(angle):
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s], [s, c]])


def _homlib_problem(seed=42, n_points=30):
    """Generate a homlib problem instance and extract ground-truth fields."""
    cfg = homlib.ProblemConfig(homlib.DistortionCase.NO_DISTORTION, 0.0, n_points)
    cfg.seed = seed
    prob = homlib.generate_problem_instance(cfg)
    x1 = np.array(prob.x1)
    x2 = np.array(prob.x2)
    A_rel = np.array(prob.A)
    ori = np.array(prob.ori)
    q = np.sqrt(np.abs(np.linalg.det(A_rel)))
    return x1, x2, A_rel, ori, q, prob


def _homlib_affine_result(seed=42, n_points=30):
    """Build an AffineMatchResult from homlib's ground-truth (ori, q).

    Individual frames are similarity transforms: affines1 = R(theta1),
    affines2 = q * R(theta2).  Their relative affine is the similarity
    approximation of A, not A itself (which also has shear).
    """
    x1, x2, A_rel, ori, q, prob = _homlib_problem(seed, n_points)
    n = len(x1)
    theta1, theta2 = ori[:, 0], ori[:, 1]
    affines1 = np.array([_rotation_matrix(a) for a in theta1])
    affines2 = np.array([qi * _rotation_matrix(a) for qi, a in zip(q, theta2)])
    return AffineMatchResult(
        points1=x1, points2=x2, scores=None,
        angles1=theta1, angles2=theta2,
        scales1=np.ones(n), scales2=q,
        affines1=affines1, affines2=affines2,
    ), A_rel, ori, q, prob


class TestAffineDecompositionUpgrader(unittest.TestCase):
    """Downgrade: AffineMatchResult -> ScaleOrientationMatchResult.

    homlib's affine2sift decomposes A into (theta1, theta2, q).  Individual
    similarity frames s*R(theta) decompose EXACTLY via atan2/det, so the
    round-trip through AffineDecompositionUpgrader must recover the ground-truth
    angles and scales to machine precision.
    """

    def test_homlib_angles_exact(self):
        """Decomposing similarity frames built from homlib's (ori, q)
        recovers the ground-truth theta1, theta2 exactly."""
        result, A_rel, ori, q, _ = _homlib_affine_result(seed=42, n_points=40)
        out = AffineDecompositionUpgrader().upgrade(result)

        self.assertIsInstance(out, ScaleOrientationMatchResult)
        self.assertNotIsInstance(out, AffineMatchResult)
        np.testing.assert_allclose(out.angles1, ori[:, 0], atol=1e-12)
        np.testing.assert_allclose(out.angles2, ori[:, 1], atol=1e-12)

    def test_homlib_scales_exact(self):
        """Decomposing similarity frames recovers scales1=1 and scales2=q exactly."""
        result, A_rel, ori, q, _ = _homlib_affine_result(seed=99, n_points=40)
        out = AffineDecompositionUpgrader().upgrade(result)

        np.testing.assert_allclose(out.scales1, 1.0, atol=1e-12)
        np.testing.assert_allclose(out.scales2, q, atol=1e-12)

    def test_homlib_solver_consistency(self):
        """The decomposed angles feed correctly into the orientation solver,
        recovering the ground-truth homography."""
        result, A_rel, ori, q, prob = _homlib_affine_result(seed=77, n_points=10)

        solver = homlib.ValtonenOrnhagICPR2026OrientationOneSided()
        poses = solver.estimate(
            result.points1.T, result.points2.T,
            ori=np.array(prob.ori).T,
        )
        self.assertGreater(len(poses), 0, "Orientation solver should find a pose")
        self.assertLess(prob.hom_error(poses[0].homography), 1e-6)

    def test_multiple_seeds(self):
        """Exact recovery holds across many random problem instances."""
        for seed in range(10):
            result, _, ori, q, _ = _homlib_affine_result(seed=seed, n_points=20)
            out = AffineDecompositionUpgrader().upgrade(result)
            np.testing.assert_allclose(out.angles1, ori[:, 0], atol=1e-12,
                                       err_msg=f"angles1 failed for seed={seed}")
            np.testing.assert_allclose(out.angles2, ori[:, 1], atol=1e-12,
                                       err_msg=f"angles2 failed for seed={seed}")
            np.testing.assert_allclose(out.scales2, q, atol=1e-12,
                                       err_msg=f"scales2 failed for seed={seed}")


class TestBarathAffineUpgrader(unittest.TestCase):
    """Upgrade: ScaleOrientationMatchResult -> AffineMatchResult (approximate).

    The Barath upgrader recovers the full relative affine from scale+orientation
    and the fundamental matrix F.  The scene is non-planar (unlike homlib's
    generator), so we construct synthetic cameras with known F.
    """

    @staticmethod
    def _make_two_view_geometry(rng, n=30):
        """Two calibrated cameras with known F and non-planar 3D points."""
        K = np.diag([500.0, 500.0, 1.0])
        R3 = np.eye(3)
        R3[:2, :2] = _rotation_matrix(0.3)
        t = np.array([0.5, 0.1, 0.2])
        tx = np.array([[0, -t[2], t[1]], [t[2], 0, -t[0]], [-t[1], t[0], 0]])
        E = tx @ R3
        K_inv = np.linalg.inv(K)
        F = K_inv.T @ E @ K_inv
        F /= np.linalg.norm(F)

        pts3d = rng.uniform(-1, 1, (n, 3))
        pts3d[:, 2] += 5.0
        P1 = K @ np.hstack([np.eye(3), np.zeros((3, 1))])
        P2 = K @ np.hstack([R3, t.reshape(3, 1)])
        x1h = (P1 @ np.hstack([pts3d, np.ones((n, 1))]).T).T
        x2h = (P2 @ np.hstack([pts3d, np.ones((n, 1))]).T).T
        x1 = x1h[:, :2] / x1h[:, 2:3]
        x2 = x2h[:, :2] / x2h[:, 2:3]
        return F, x1, x2

    def test_affines1_exact(self):
        """The output affines1 should be exactly s1*R(a1) (similarity by construction)."""
        rng = np.random.default_rng(789)
        F, x1, x2 = self._make_two_view_geometry(rng, n=40)
        n = len(x1)
        # Use homlib-like angles and scales
        x1h, _, _, ori, q, _ = _homlib_problem(seed=42, n_points=n)
        angles1, angles2 = ori[:, 0], ori[:, 1]
        scales1, scales2 = np.ones(n), q

        inp = ScaleOrientationMatchResult(
            points1=x1, points2=x2, scores=None,
            angles1=angles1, angles2=angles2,
            scales1=scales1, scales2=scales2,
        )
        out = BarathAffineUpgrader().upgrade(inp, F=F)

        self.assertIsInstance(out, AffineMatchResult)
        for i in range(n):
            expected = scales1[i] * _rotation_matrix(angles1[i])
            np.testing.assert_allclose(out.affines1[i], expected, atol=1e-10)

    def test_relative_affine_det_positive(self):
        """Reconstructed relative affines have positive determinant and finite entries."""
        rng = np.random.default_rng(321)
        F, x1, x2 = self._make_two_view_geometry(rng, n=30)
        n = len(x1)
        _, _, _, ori, q, _ = _homlib_problem(seed=99, n_points=n)

        inp = ScaleOrientationMatchResult(
            points1=x1, points2=x2, scores=None,
            angles1=ori[:, 0], angles2=ori[:, 1],
            scales1=np.ones(n), scales2=q,
        )
        out = BarathAffineUpgrader().upgrade(inp, F=F)
        A_rel = out.relative_affines()

        self.assertTrue(np.all(np.isfinite(A_rel)))
        dets = np.linalg.det(A_rel)
        self.assertTrue(np.all(dets > 0))

    def test_downgrade_round_trip(self):
        """Upgrade then downgrade: decomposing affines1 (= s1*R1) recovers the
        original angles and scales exactly."""
        rng = np.random.default_rng(555)
        F, x1, x2 = self._make_two_view_geometry(rng, n=25)
        n = len(x1)
        _, _, _, ori, q, _ = _homlib_problem(seed=77, n_points=n)
        angles1, angles2 = ori[:, 0], ori[:, 1]
        scales1, scales2 = np.ones(n), q

        inp = ScaleOrientationMatchResult(
            points1=x1, points2=x2, scores=None,
            angles1=angles1, angles2=angles2,
            scales1=scales1, scales2=scales2,
        )
        upgraded = BarathAffineUpgrader().upgrade(inp, F=F)
        downgraded = AffineDecompositionUpgrader().upgrade(upgraded)

        # affines1 = s1 * R1, so decomposition recovers angles1 and scales1 exactly
        np.testing.assert_allclose(downgraded.angles1, angles1, atol=1e-10)
        np.testing.assert_allclose(downgraded.scales1, scales1, atol=1e-10)

    def test_scale_ratio_preserved(self):
        """The determinant of the reconstructed relative affine should
        approximately equal (q2/q1)^2, the squared scale ratio."""
        rng = np.random.default_rng(111)
        F, x1, x2 = self._make_two_view_geometry(rng, n=30)
        n = len(x1)
        _, _, _, ori, q, _ = _homlib_problem(seed=55, n_points=n)

        inp = ScaleOrientationMatchResult(
            points1=x1, points2=x2, scores=None,
            angles1=ori[:, 0], angles2=ori[:, 1],
            scales1=np.ones(n), scales2=q,
        )
        out = BarathAffineUpgrader().upgrade(inp, F=F)
        A_rel = out.relative_affines()

        dets = np.linalg.det(A_rel)
        expected_dets = (q / 1.0) ** 2
        np.testing.assert_allclose(dets, expected_dets, rtol=0.3,
                                   err_msg="det(A_rel) should approximate (q2/q1)^2")


if __name__ == '__main__':
    unittest.main()
