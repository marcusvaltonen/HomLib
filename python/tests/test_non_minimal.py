"""Tests for NonMinimalSettings validation and edge cases."""
import unittest

import homlib


class ParseNonMinimalStrategyTestCase(unittest.TestCase):
    """Tests for _parse_non_minimal_strategy (exposed through constructors)."""

    def test_string_sample(self):
        lo = homlib.LOMSAC(homlib.DLT(), non_minimal_strategy="sample",
                           inlier_threshold=1.0)
        self.assertEqual(lo.non_minimal_strategy,
                         homlib.NonMinimalStrategy.Sample)

    def test_string_overdetermined(self):
        lo = homlib.LOMSAC(homlib.FitzgibbonCVPR2001OneSided(),
                           non_minimal_strategy="overdetermined",
                           inlier_threshold=1.0)
        self.assertEqual(lo.non_minimal_strategy,
                         homlib.NonMinimalStrategy.Overdetermined)

    def test_string_case_insensitive(self):
        lo = homlib.LOMSAC(homlib.FitzgibbonCVPR2001OneSided(),
                           non_minimal_strategy="Overdetermined",
                           inlier_threshold=1.0)
        self.assertEqual(lo.non_minimal_strategy,
                         homlib.NonMinimalStrategy.Overdetermined)

    def test_enum_value(self):
        lo = homlib.LOMSAC(homlib.DLT(),
                           non_minimal_strategy=homlib.NonMinimalStrategy.Sample,
                           inlier_threshold=1.0)
        self.assertEqual(lo.non_minimal_strategy,
                         homlib.NonMinimalStrategy.Sample)

    def test_invalid_string_raises(self):
        with self.assertRaises(ValueError):
            homlib.LOMSAC(homlib.DLT(), non_minimal_strategy="bogus",
                          inlier_threshold=1.0)

    def test_invalid_type_raises(self):
        with self.assertRaises(TypeError):
            homlib.LOMSAC(homlib.DLT(), non_minimal_strategy=42,
                          inlier_threshold=1.0)

    def test_none_defaults_to_sample(self):
        lo = homlib.LOMSAC(homlib.DLT(), inlier_threshold=1.0)
        self.assertEqual(lo.non_minimal_strategy,
                         homlib.NonMinimalStrategy.Sample)


class CheckOverdeterminedTestCase(unittest.TestCase):
    """Tests for _check_overdetermined (solver must support OD)."""

    def test_od_with_non_od_solver_raises(self):
        with self.assertRaises(ValueError, msg="Wadenback does not support OD"):
            homlib.LOMSAC(homlib.Wadenback3DV2026OneSided(),
                          non_minimal_strategy="overdetermined",
                          inlier_threshold=1.0)

    def test_od_with_od_solver_ok(self):
        lo = homlib.LOMSAC(homlib.FitzgibbonCVPR2001OneSided(),
                           non_minimal_strategy="overdetermined",
                           inlier_threshold=1.0)
        self.assertEqual(lo.non_minimal_strategy,
                         homlib.NonMinimalStrategy.Overdetermined)

    def test_od_with_dlt_ok(self):
        lo = homlib.LOMSAC(homlib.DLT(),
                           non_minimal_strategy="overdetermined",
                           inlier_threshold=1.0)
        self.assertEqual(lo.non_minimal_strategy,
                         homlib.NonMinimalStrategy.Overdetermined)


class ResolveNonMinimalSolverTestCase(unittest.TestCase):
    """Tests for _resolve_non_minimal_solver."""

    def test_solver_object(self):
        lo = homlib.LOMSAC(
            homlib.Wadenback3DV2026OneSided(),
            non_minimal_solver=homlib.FitzgibbonCVPR2001OneSided(),
            inlier_threshold=1.0)
        self.assertIsNotNone(lo.non_minimal_solver)

    def test_solver_id(self):
        lo = homlib.LOMSAC(
            homlib.Wadenback3DV2026OneSided(),
            non_minimal_solver=homlib.SolverId.FitzgibbonCVPR2001OneSided,
            inlier_threshold=1.0)
        self.assertIsNotNone(lo.non_minimal_solver)

    def test_string_name(self):
        lo = homlib.LOMSAC(
            homlib.Wadenback3DV2026OneSided(),
            non_minimal_solver="FitzgibbonCVPR2001OneSided",
            inlier_threshold=1.0)
        self.assertIsNotNone(lo.non_minimal_solver)

    def test_auto_sets_overdetermined(self):
        lo = homlib.LOMSAC(
            homlib.Wadenback3DV2026OneSided(),
            non_minimal_solver=homlib.FitzgibbonCVPR2001OneSided(),
            inlier_threshold=1.0)
        self.assertEqual(lo.non_minimal_strategy,
                         homlib.NonMinimalStrategy.Overdetermined)

    def test_mismatched_distortion_case_raises(self):
        with self.assertRaises(ValueError, msg="distortion_case mismatch"):
            homlib.LOMSAC(
                homlib.Wadenback3DV2026TwoSided(),
                non_minimal_solver=homlib.FitzgibbonCVPR2001OneSided(),
                inlier_threshold=1.0)

    def test_affine_nm_solver_raises(self):
        with self.assertRaises(ValueError, msg="affine solver not point-based"):
            homlib.LOMSAC(
                homlib.Wadenback3DV2026OneSided(),
                non_minimal_solver=homlib.BarathVISAPP2016Affine(),
                inlier_threshold=1.0)

    def test_non_od_nm_solver_raises(self):
        with self.assertRaises(ValueError, msg="nm solver must support OD"):
            homlib.LOMSAC(
                homlib.Wadenback3DV2026OneSided(),
                non_minimal_solver=homlib.Wadenback3DV2026OneSided(),
                inlier_threshold=1.0)

    def test_wrong_type_raises(self):
        with self.assertRaises(TypeError):
            homlib.LOMSAC(
                homlib.DLT(),
                non_minimal_solver=42,
                inlier_threshold=1.0)


class CheckNonMinimalParamsTestCase(unittest.TestCase):
    """Tests for _check_non_minimal_params."""

    def test_nm_solver_with_sample_strategy_raises(self):
        with self.assertRaises(ValueError, msg="nm_solver requires OD strategy"):
            homlib.LOMSAC(
                homlib.Wadenback3DV2026OneSided(),
                non_minimal_solver=homlib.FitzgibbonCVPR2001OneSided(),
                non_minimal_strategy="sample",
                inlier_threshold=1.0)

    def test_nm_solver_with_od_strategy_ok(self):
        lo = homlib.LOMSAC(
            homlib.Wadenback3DV2026OneSided(),
            non_minimal_solver=homlib.FitzgibbonCVPR2001OneSided(),
            non_minimal_strategy="overdetermined",
            inlier_threshold=1.0)
        self.assertIsNotNone(lo.non_minimal_solver)

    def test_od_checks_nm_solver_not_primary(self):
        with self.assertRaises(ValueError):
            homlib.LOMSAC(
                homlib.Wadenback3DV2026OneSided(),
                non_minimal_solver=homlib.Wadenback3DV2026OneSided(),
                non_minimal_strategy="overdetermined",
                inlier_threshold=1.0)


@unittest.skipUnless(homlib._HAS_GCRANSAC, "GC-RANSAC not built")
class GCRANSACValidationTestCase(unittest.TestCase):
    """Same validation rules apply to GC-RANSAC."""

    def test_od_with_non_od_solver_raises(self):
        with self.assertRaises(ValueError):
            homlib.GCRANSAC(homlib.Wadenback3DV2026OneSided(),
                            non_minimal_strategy="overdetermined",
                            threshold=1.0)

    def test_nm_solver_with_sample_raises(self):
        with self.assertRaises(ValueError):
            homlib.GCRANSAC(
                homlib.Wadenback3DV2026OneSided(),
                non_minimal_solver=homlib.FitzgibbonCVPR2001OneSided(),
                non_minimal_strategy="sample",
                threshold=1.0)

    def test_nm_solver_auto_sets_od(self):
        gc = homlib.GCRANSAC(
            homlib.Wadenback3DV2026OneSided(),
            non_minimal_solver=homlib.FitzgibbonCVPR2001OneSided(),
            threshold=1.0)
        self.assertEqual(gc.non_minimal_strategy,
                         homlib.NonMinimalStrategy.Overdetermined)


@unittest.skipUnless(homlib._HAS_MAGSAC, "MAGSAC not built")
class MAGSACValidationTestCase(unittest.TestCase):
    """Same validation rules apply to MAGSAC."""

    def test_od_with_non_od_solver_raises(self):
        with self.assertRaises(ValueError):
            homlib.MAGSAC(homlib.Wadenback3DV2026OneSided(),
                          non_minimal_strategy="overdetermined",
                          sigma_max=1.0)

    def test_nm_solver_auto_sets_od(self):
        mg = homlib.MAGSAC(
            homlib.Wadenback3DV2026OneSided(),
            non_minimal_solver=homlib.FitzgibbonCVPR2001OneSided(),
            sigma_max=1.0)
        self.assertEqual(mg.non_minimal_strategy,
                         homlib.NonMinimalStrategy.Overdetermined)


@unittest.skipUnless(homlib._HAS_SUPERANSAC, "SupeRANSAC not built")
class SupeRANSACValidationTestCase(unittest.TestCase):
    """Same validation rules apply to SupeRANSAC."""

    def test_od_with_non_od_solver_raises(self):
        with self.assertRaises(ValueError):
            homlib.SupeRANSAC(homlib.Wadenback3DV2026OneSided(),
                              non_minimal_strategy="overdetermined",
                              inlier_threshold=1.0)

    def test_nm_solver_auto_sets_od(self):
        sa = homlib.SupeRANSAC(
            homlib.Wadenback3DV2026OneSided(),
            non_minimal_solver=homlib.FitzgibbonCVPR2001OneSided(),
            inlier_threshold=1.0)
        self.assertEqual(sa.non_minimal_strategy,
                         homlib.NonMinimalStrategy.Overdetermined)


class LONonMinimalStrategyPropertyTestCase(unittest.TestCase):
    """Test the non_minimal_strategy property setter."""

    def test_setter_accepts_string(self):
        lo = homlib.LOMSAC(homlib.FitzgibbonCVPR2001OneSided(),
                           inlier_threshold=1.0)
        lo.non_minimal_strategy = "overdetermined"
        self.assertEqual(lo.non_minimal_strategy,
                         homlib.NonMinimalStrategy.Overdetermined)

    def test_setter_rejects_incompatible(self):
        lo = homlib.LOMSAC(homlib.DLT(),
                           non_minimal_solver=homlib.DLT(),
                           inlier_threshold=1.0)
        with self.assertRaises(ValueError):
            lo.non_minimal_strategy = "sample"
