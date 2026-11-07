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


if __name__ == "__main__":
    unittest.main()
