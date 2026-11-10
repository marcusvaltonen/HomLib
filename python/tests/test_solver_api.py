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


class DistortionCaseMismatchTestCase(unittest.TestCase):
    """Distortion-case mismatch raises ValueError at construction."""

    MISMATCHES = [
        ("DLT", "NakanoICPR2025OneSided"),
        ("DLT", "FitzgibbonCVPR2001TwoSidedEqual"),
        ("FitzgibbonCVPR2001OneSided", "DLT"),
        ("KukelovaCVPR2015TwoSided", "NakanoICPR2025OneSided"),
        ("Wadenback3DV2026TwoSidedEqual", "NakanoICPR2025OneSidedRight"),
    ]

    def test_lomsac_mismatch(self):
        for minimal, non_minimal in self.MISMATCHES:
            with self.subTest(minimal=minimal, non_minimal=non_minimal):
                with self.assertRaises(ValueError) as ctx:
                    homlib.LOMSAC(
                        minimal,
                        non_minimal_strategy="overdetermined",
                        non_minimal_solver=non_minimal,
                        inlier_threshold=3.0,
                    )
                self.assertIn("distortion cases must match", str(ctx.exception))

    def test_superansac_mismatch(self):
        with self.assertRaises(ValueError):
            homlib.SupeRANSAC(
                "DLT",
                non_minimal_strategy="overdetermined",
                non_minimal_solver="NakanoICPR2025OneSided",
            )

    def test_gcransac_mismatch(self):
        with self.assertRaises(ValueError):
            homlib.GCRANSAC(
                "DLT",
                non_minimal_strategy="overdetermined",
                non_minimal_solver="NakanoICPR2025OneSided",
            )

    def test_magsac_mismatch(self):
        with self.assertRaises(ValueError):
            homlib.MAGSAC(
                "DLT",
                non_minimal_strategy="overdetermined",
                non_minimal_solver="NakanoICPR2025OneSided",
            )


class DistortionCaseMatchTestCase(unittest.TestCase):
    """Compatible non-minimal solver combinations succeed."""

    COMPATIBLE = [
        ("CaiPAMI2025", "DLT"),
        ("Wadenback3DV2026OneSided", "NakanoICPR2025OneSided"),
        ("Wadenback3DV2026OneSided", "FitzgibbonCVPR2001OneSided"),
        ("Wadenback3DV2026OneSidedRight", "NakanoICPR2025OneSidedRight"),
        ("KukelovaCVPR2015TwoSidedEqual", "FitzgibbonCVPR2001TwoSidedEqual"),
    ]

    def test_lomsac_compatible(self):
        for minimal, non_minimal in self.COMPATIBLE:
            with self.subTest(minimal=minimal, non_minimal=non_minimal):
                est = homlib.LOMSAC(
                    minimal,
                    non_minimal_strategy="overdetermined",
                    non_minimal_solver=non_minimal,
                    inlier_threshold=3.0,
                )
                self.assertIsNotNone(est)

    def test_non_minimal_solver_as_enum(self):
        est = homlib.LOMSAC(
            homlib.SolverId.Wadenback3DV2026OneSided,
            non_minimal_strategy="overdetermined",
            non_minimal_solver=homlib.SolverId.NakanoICPR2025OneSided,
            inlier_threshold=3.0,
        )
        self.assertIsNotNone(est)

    def test_non_minimal_solver_as_solver_instance(self):
        est = homlib.LOMSAC(
            "CaiPAMI2025",
            non_minimal_strategy="overdetermined",
            non_minimal_solver=homlib.DLT(use_qr=True),
            inlier_threshold=3.0,
        )
        self.assertIsNotNone(est)

    def test_non_minimal_solver_instance_preserves_type(self):
        nm = homlib.DLT(use_qr=True)
        est = homlib.LOMSAC(
            "CaiPAMI2025",
            non_minimal_strategy="overdetermined",
            non_minimal_solver=nm,
            inlier_threshold=3.0,
        )
        self.assertTrue(est._non_minimal_solver.options["use_qr"])


class NonMinimalStrategyValidationTestCase(unittest.TestCase):
    """non_minimal_solver without overdetermined strategy raises ValueError."""

    def test_lomsac_default_strategy(self):
        with self.assertRaises(ValueError) as ctx:
            homlib.LOMSAC(
                "Wadenback3DV2026OneSided",
                non_minimal_solver="NakanoICPR2025OneSided",
                inlier_threshold=3.0,
            )
        self.assertIn("overdetermined", str(ctx.exception))

    def test_lomsac_explicit_sample(self):
        with self.assertRaises(ValueError):
            homlib.LOMSAC(
                "Wadenback3DV2026OneSided",
                non_minimal_strategy="sample",
                non_minimal_solver="NakanoICPR2025OneSided",
                inlier_threshold=3.0,
            )

    def test_superansac_default_strategy(self):
        with self.assertRaises(ValueError):
            homlib.SupeRANSAC(
                "Wadenback3DV2026OneSided",
                non_minimal_solver="NakanoICPR2025OneSided",
            )

    def test_gcransac_default_strategy(self):
        with self.assertRaises(ValueError):
            homlib.GCRANSAC(
                "Wadenback3DV2026OneSided",
                non_minimal_solver="NakanoICPR2025OneSided",
            )

    def test_magsac_default_strategy(self):
        with self.assertRaises(ValueError):
            homlib.MAGSAC(
                "Wadenback3DV2026OneSided",
                non_minimal_solver="NakanoICPR2025OneSided",
            )

    def test_lomsac_overdetermined_succeeds(self):
        est = homlib.LOMSAC(
            "Wadenback3DV2026OneSided",
            non_minimal_strategy="overdetermined",
            non_minimal_solver="NakanoICPR2025OneSided",
            inlier_threshold=3.0,
        )
        self.assertIsNotNone(est)

    def test_no_solver_no_error(self):
        """strategy=sample without non_minimal_solver is fine."""
        est = homlib.LOMSAC("DLT", non_minimal_strategy="sample", inlier_threshold=3.0)
        self.assertIsNotNone(est)

    def test_overdetermined_without_solver_is_fine(self):
        """strategy=overdetermined without non_minimal_solver falls back to sample in C++."""
        est = homlib.LOMSAC("DLT", non_minimal_strategy="overdetermined", inlier_threshold=3.0)
        self.assertIsNotNone(est)


class StringInBackendsTestCase(unittest.TestCase):
    """All backends accept string solver names."""

    def test_lomsac(self):
        est = homlib.LOMSAC("DLT", inlier_threshold=3.0)
        self.assertEqual(est.solver.id, homlib.SolverId.DLT)

    def test_superansac(self):
        est = homlib.SupeRANSAC("DLT")
        self.assertEqual(est.solver.id, homlib.SolverId.DLT)

    def test_gcransac(self):
        est = homlib.GCRANSAC("DLT")
        self.assertEqual(est.solver.id, homlib.SolverId.DLT)

    def test_magsac(self):
        est = homlib.MAGSAC("DLT")
        self.assertEqual(est.solver.id, homlib.SolverId.DLT)


class SupportsOverdeterminedTestCase(unittest.TestCase):
    """Solver.supports_overdetermined property and validation."""

    OVERDETERMINED = {
        sid.name for sid in homlib.SolverId
        if homlib.Solver(sid).supports_overdetermined
    }

    def test_overdetermined_set_size(self):
        self.assertEqual(len(self.OVERDETERMINED), 6)

    def test_supports_overdetermined_true(self):
        for name in self.OVERDETERMINED:
            with self.subTest(name=name):
                self.assertTrue(homlib.Solver(name).supports_overdetermined)

    def test_supports_overdetermined_false(self):
        for sid in homlib.SolverId:
            if sid.name not in self.OVERDETERMINED:
                with self.subTest(name=sid.name):
                    self.assertFalse(homlib.Solver(sid).supports_overdetermined)

    def test_all_non_overdetermined_solvers_rejected(self):
        """Every non-overdetermined point solver is rejected as non-minimal."""
        point_solvers = homlib.solver_ids(homlib.FeatureType.Points)
        non_od = [sid for sid in point_solvers if sid.name not in self.OVERDETERMINED]
        self.assertTrue(len(non_od) > 0)
        for sid in non_od:
            dc = homlib.Solver(sid).distortion_case
            minimal = next(
                s for s in point_solvers
                if homlib.Solver(s).distortion_case == dc and s != sid
            )
            with self.subTest(solver=sid.name):
                with self.assertRaises(ValueError) as ctx:
                    homlib.LOMSAC(
                        minimal,
                        non_minimal_strategy="overdetermined",
                        non_minimal_solver=sid,
                        inlier_threshold=3.0,
                    )
                self.assertIn("does not support overdetermined", str(ctx.exception))

    def test_error_lists_compatible_solvers(self):
        with self.assertRaises(ValueError) as ctx:
            homlib.LOMSAC(
                "Wadenback3DV2026OneSided",
                non_minimal_strategy="overdetermined",
                non_minimal_solver="Wadenback3DV2026OneSided",
                inlier_threshold=3.0,
            )
        msg = str(ctx.exception)
        self.assertIn("does not support overdetermined", msg)
        self.assertIn("FitzgibbonCVPR2001OneSided", msg)
        self.assertIn("NakanoICPR2025OneSided", msg)

    def test_no_distortion_compatible_list(self):
        with self.assertRaises(ValueError) as ctx:
            homlib.LOMSAC(
                "DLT",
                non_minimal_strategy="overdetermined",
                non_minimal_solver="CaiPAMI2025",
                inlier_threshold=3.0,
            )
        msg = str(ctx.exception)
        self.assertIn("DLT", msg)

    def test_sample_strategy_no_overdetermined_check(self):
        """Sample strategy without non_minimal_solver needs no validation."""
        est = homlib.LOMSAC("DLT", non_minimal_strategy="sample", inlier_threshold=3.0)
        self.assertIsNotNone(est)

    def test_readonly(self):
        s = homlib.Solver("DLT")
        with self.assertRaises(AttributeError):
            s.supports_overdetermined = False


if __name__ == "__main__":
    unittest.main()
