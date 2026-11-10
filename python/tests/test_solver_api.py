import unittest

import homlib


class DistortionCasePropertyTestCase(unittest.TestCase):
    """Solver.distortion_case: readonly property from C++."""

    EXPECTED = {
        homlib.SolverId.DLT: homlib.DistortionCase.NO_DISTORTION,
        homlib.SolverId.CaiPAMI2025: homlib.DistortionCase.NO_DISTORTION,
        homlib.SolverId.GuoPRL2019: homlib.DistortionCase.NO_DISTORTION,
        homlib.SolverId.FitzgibbonCVPR2001OneSided: homlib.DistortionCase.ONE_SIDED_LEFT,
        homlib.SolverId.NakanoICPR2025OneSided: homlib.DistortionCase.ONE_SIDED_LEFT,
        homlib.SolverId.Wadenback3DV2026OneSided: homlib.DistortionCase.ONE_SIDED_LEFT,
        homlib.SolverId.FitzgibbonCVPR2001OneSidedRight: homlib.DistortionCase.ONE_SIDED_RIGHT,
        homlib.SolverId.NakanoICPR2025OneSidedRight: homlib.DistortionCase.ONE_SIDED_RIGHT,
        homlib.SolverId.Wadenback3DV2026OneSidedRight: homlib.DistortionCase.ONE_SIDED_RIGHT,
        homlib.SolverId.FitzgibbonCVPR2001TwoSidedEqual: homlib.DistortionCase.TWO_SIDED_EQUAL,
        homlib.SolverId.KukelovaCVPR2015TwoSidedEqual: homlib.DistortionCase.TWO_SIDED_EQUAL,
        homlib.SolverId.KukelovaCVPR2015TwoSidedEqual6Pt: homlib.DistortionCase.TWO_SIDED_EQUAL,
        homlib.SolverId.Wadenback3DV2026TwoSidedEqual: homlib.DistortionCase.TWO_SIDED_EQUAL,
        homlib.SolverId.KukelovaCVPR2015TwoSided: homlib.DistortionCase.TWO_SIDED,
        homlib.SolverId.KukelovaCVPR2015TwoSided6Pt: homlib.DistortionCase.TWO_SIDED,
        homlib.SolverId.Wadenback3DV2026TwoSided: homlib.DistortionCase.TWO_SIDED,
        homlib.SolverId.BarathVISAPP2016Affine: homlib.DistortionCase.NO_DISTORTION,
        homlib.SolverId.ValtonenOrnhagICPR2026AffineOneSided: homlib.DistortionCase.ONE_SIDED_LEFT,
        homlib.SolverId.ValtonenOrnhagICPR2026AffineOneSidedRight: homlib.DistortionCase.ONE_SIDED_RIGHT,
        homlib.SolverId.ValtonenOrnhagICPR2026OrientationOneSided: homlib.DistortionCase.ONE_SIDED_LEFT,
        homlib.SolverId.ValtonenOrnhagICPR2026OrientationOneSidedRight: homlib.DistortionCase.ONE_SIDED_RIGHT,
    }

    def test_representative_solvers(self):
        for solver_id, expected_dc in self.EXPECTED.items():
            with self.subTest(solver=solver_id.name):
                self.assertEqual(homlib.Solver(solver_id).distortion_case, expected_dc)

    def test_all_solvers_covered(self):
        """Every SolverId is in the EXPECTED dict — no solver left untested."""
        self.assertEqual(set(self.EXPECTED.keys()), set(homlib.SolverId))

    def test_all_solvers_have_distortion_case(self):
        dc_values = set(homlib.DistortionCase)
        for solver_id in homlib.SolverId:
            with self.subTest(solver=solver_id.name):
                dc = homlib.Solver(solver_id).distortion_case
                self.assertIn(dc, dc_values)

    def test_distortion_case_readonly(self):
        s = homlib.Solver(homlib.SolverId.DLT)
        with self.assertRaises(AttributeError):
            s.distortion_case = homlib.DistortionCase.TWO_SIDED


if __name__ == "__main__":
    unittest.main()
