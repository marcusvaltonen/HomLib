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


class SolverStringConstructionTestCase(unittest.TestCase):
    """Solver("name") string convenience."""

    def test_from_string(self):
        s = homlib.Solver("DLT")
        self.assertEqual(s.id, homlib.SolverId.DLT)

    def test_from_string_long_name(self):
        s = homlib.Solver("NakanoICPR2025OneSided")
        self.assertEqual(s.id, homlib.SolverId.NakanoICPR2025OneSided)

    def test_all_solvers_from_string(self):
        for solver_id in homlib.SolverId:
            with self.subTest(name=solver_id.name):
                s = homlib.Solver(solver_id.name)
                self.assertEqual(s.id, solver_id)

    def test_from_string_invalid(self):
        with self.assertRaises(ValueError):
            homlib.Solver("NotASolver")

    def test_from_solver_id(self):
        s = homlib.Solver(homlib.SolverId.DLT)
        self.assertEqual(s.id, homlib.SolverId.DLT)

    def test_string_preserves_kwargs(self):
        s = homlib.Solver("DLT", use_qr=True)
        self.assertEqual(s.id, homlib.SolverId.DLT)
        self.assertTrue(s.options["use_qr"])

    def test_string_in_named_solver(self):
        """Named solver classes still work (they use SolverId internally)."""
        s = homlib.NakanoICPR2025OneSided(extra_check=True)
        self.assertEqual(s.id, homlib.SolverId.NakanoICPR2025OneSided)
        self.assertTrue(s.options["extra_check"])


class StringInBackendsTestCase(unittest.TestCase):
    """All backends accept string solver names."""

    def test_lomsac(self):
        est = homlib.LOMSAC("DLT", inlier_threshold=3.0)
        self.assertEqual(est.solver.id, homlib.SolverId.DLT)

    @unittest.skipUnless(homlib._HAS_SUPERANSAC, "SupeRANSAC not built")
    def test_superansac(self):
        est = homlib.SupeRANSAC("DLT")
        self.assertEqual(est.solver.id, homlib.SolverId.DLT)

    @unittest.skipUnless(homlib._HAS_GCRANSAC, "GC-RANSAC not built")
    def test_gcransac(self):
        est = homlib.GCRANSAC("DLT")
        self.assertEqual(est.solver.id, homlib.SolverId.DLT)

    @unittest.skipUnless(homlib._HAS_MAGSAC, "MAGSAC not built")
    def test_magsac(self):
        est = homlib.MAGSAC("DLT")
        self.assertEqual(est.solver.id, homlib.SolverId.DLT)


if __name__ == "__main__":
    unittest.main()
