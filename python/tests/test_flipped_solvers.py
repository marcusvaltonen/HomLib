import unittest

import numpy as np
import homlib


class FlippedSolverConsistencyTestCase(unittest.TestCase):
    """Flipped (left/default) vs native (right) solvers:
    inverted H and swapped distortion parameters."""

    def setUp(self):
        np.random.seed(42)
        self.p1 = np.random.randn(2, 5) * 0.5
        self.p2 = np.random.randn(2, 5) * 0.5

    def _verify_consistency(self, right_cls, left_cls):
        right_results = right_cls().estimate(self.p1, self.p2)
        left_results = left_cls().estimate(self.p2, self.p1)

        self.assertEqual(len(right_results), len(left_results))
        self.assertGreater(len(right_results), 0)

        for r, l in zip(right_results, left_results):
            self.assertAlmostEqual(
                l.distortion_parameter,
                r.distortion_parameter2, places=10,
            )
            self.assertAlmostEqual(
                l.distortion_parameter2,
                r.distortion_parameter, places=10,
            )
            H_inv = np.linalg.inv(r.homography)
            np.testing.assert_allclose(
                l.homography, H_inv, atol=1e-10,
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
    """Flipped affine solver: left vs native right."""

    def test_affine(self):
        x = np.array([
            [-0.304372059145248, -0.308082306270758],
            [-0.174125682172663, 0.100342577802908],
        ])
        y = np.array([
            [0.596062339231747, 0.284590553127076],
            [0.449851750693698, 0.20125992440591],
        ])
        A1 = np.array([
            [-0.742737063351745, -1.25720051435581],
            [1.29931162730836, -0.953395533218031],
        ])
        A2 = np.array([
            [-0.763314941736403, -1.05311577674003],
            [1.13911319830585, -0.839438880618645],
        ])
        A = np.column_stack([A1.flatten(), A2.flatten()])
        A_inv = np.column_stack([
            np.linalg.inv(A1).flatten(),
            np.linalg.inv(A2).flatten(),
        ])

        right = homlib.ValtonenOrnhagICPR2026AffineOneSidedRight()
        rr = right.estimate(x, y, A=A)
        left = homlib.ValtonenOrnhagICPR2026AffineOneSided()
        lr = left.estimate(y, x, A=A_inv)

        self.assertEqual(len(rr), len(lr))
        self.assertGreater(len(rr), 0)
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
    """Flipped orientation solver: left vs native right."""

    def test_orientation(self):
        x = np.array([
            [0.611733207361777, 0.506500688761625,
             0.618492431042353, 0.446207436531052],
            [-0.100815813483642, 0.22654451870654,
             0.490634389234033, -0.221543454776215],
        ])
        y = np.array([
            [0.152378076324475, 0.35249479911981,
             0.576429540272504, -0.0960825885072252],
            [0.471173894309057, 0.604057206211103,
             0.546190419126011, 0.57183541266364],
        ])
        ori = np.array([
            [2.82616873462069, 2.76670870212448,
             3.03912230291483, 2.54398570702913],
            [2.11084467936878, 1.5791745485192,
             1.40242855410085, 2.11604504644583],
        ])
        ori_swapped = np.array([ori[1], ori[0]])

        right = homlib.ValtonenOrnhagICPR2026OrientationOneSidedRight()
        rr = right.estimate(x, y, ori=ori)
        left = homlib.ValtonenOrnhagICPR2026OrientationOneSided()
        lr = left.estimate(y, x, ori=ori_swapped)

        self.assertEqual(len(rr), len(lr))
        self.assertGreater(len(rr), 0)
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
