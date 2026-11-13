"""Tests for NonMinimalSettings validation and edge cases."""
import unittest

import numpy as np
import homlib


class ParseNonMinimalStrategyTestCase(unittest.TestCase):
    """Tests for strategy parsing via non_minimal_settings dict."""

    def test_string_sample(self):
        lo = homlib.LOMSAC(homlib.DLT(),
                           non_minimal_settings={"strategy": "sample"},
                           inlier_threshold=1.0)
        self.assertEqual(lo.non_minimal_strategy,
                         homlib.NonMinimalStrategy.Sample)

    def test_string_overdetermined(self):
        lo = homlib.LOMSAC(homlib.FitzgibbonCVPR2001OneSided(),
                           non_minimal_settings={"strategy": "overdetermined"},
                           inlier_threshold=1.0)
        self.assertEqual(lo.non_minimal_strategy,
                         homlib.NonMinimalStrategy.Overdetermined)

    def test_string_case_insensitive(self):
        lo = homlib.LOMSAC(homlib.FitzgibbonCVPR2001OneSided(),
                           non_minimal_settings={"strategy": "Overdetermined"},
                           inlier_threshold=1.0)
        self.assertEqual(lo.non_minimal_strategy,
                         homlib.NonMinimalStrategy.Overdetermined)

    def test_enum_value(self):
        lo = homlib.LOMSAC(homlib.DLT(),
                           non_minimal_settings={"strategy": homlib.NonMinimalStrategy.Sample},
                           inlier_threshold=1.0)
        self.assertEqual(lo.non_minimal_strategy,
                         homlib.NonMinimalStrategy.Sample)

    def test_invalid_string_raises(self):
        with self.assertRaises(ValueError):
            homlib.LOMSAC(homlib.DLT(),
                          non_minimal_settings={"strategy": "bogus"},
                          inlier_threshold=1.0)

    def test_invalid_type_raises(self):
        with self.assertRaises(TypeError):
            homlib.LOMSAC(homlib.DLT(),
                          non_minimal_settings={"strategy": 42},
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
                          non_minimal_settings={"strategy": "overdetermined"},
                          inlier_threshold=1.0)

    def test_od_with_od_solver_ok(self):
        lo = homlib.LOMSAC(homlib.FitzgibbonCVPR2001OneSided(),
                           non_minimal_settings={"strategy": "overdetermined"},
                           inlier_threshold=1.0)
        self.assertEqual(lo.non_minimal_strategy,
                         homlib.NonMinimalStrategy.Overdetermined)

    def test_od_with_dlt_ok(self):
        lo = homlib.LOMSAC(homlib.DLT(),
                           non_minimal_settings={"strategy": "overdetermined"},
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

    def test_affine_nm_for_point_primary_raises(self):
        with self.assertRaises(ValueError, msg="affine nm for point primary"):
            homlib.LOMSAC(
                homlib.Wadenback3DV2026OneSided(),
                non_minimal_solver=homlib.BarathVISAPP2016Affine(),
                inlier_threshold=1.0)

    def test_orientation_nm_for_affine_primary_raises(self):
        with self.assertRaises(ValueError, msg="cross-type: ori nm for affine"):
            homlib.LOMSAC(
                homlib.BarathVISAPP2016Affine(),
                non_minimal_solver=homlib.ValtonenOrnhagICPR2026OrientationOneSidedRight(),
                inlier_threshold=1.0)

    def test_affine_nm_for_orientation_primary_raises(self):
        with self.assertRaises(ValueError, msg="cross-type: affine nm for ori"):
            homlib.LOMSAC(
                homlib.ValtonenOrnhagICPR2026OrientationOneSidedRight(),
                non_minimal_solver=homlib.ValtonenOrnhagICPR2026AffineOneSidedRight(),
                inlier_threshold=1.0)

    def test_affine_nm_for_affine_primary_ok(self):
        lo = homlib.LOMSAC(
            homlib.BarathVISAPP2016Affine(),
            non_minimal_solver=homlib.BarathVISAPP2016Affine(),
            inlier_threshold=1.0)
        self.assertIsNotNone(lo.non_minimal_solver)

    def test_orientation_nm_for_orientation_primary_ok(self):
        lo = homlib.LOMSAC(
            homlib.ValtonenOrnhagICPR2026OrientationOneSidedRight(),
            non_minimal_solver=homlib.ValtonenOrnhagICPR2026OrientationOneSidedRight(),
            inlier_threshold=1.0)
        self.assertIsNotNone(lo.non_minimal_solver)

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
                non_minimal_settings={"strategy": "sample"},
                inlier_threshold=1.0)

    def test_nm_solver_with_od_strategy_ok(self):
        lo = homlib.LOMSAC(
            homlib.Wadenback3DV2026OneSided(),
            non_minimal_solver=homlib.FitzgibbonCVPR2001OneSided(),
            non_minimal_settings={"strategy": "overdetermined"},
            inlier_threshold=1.0)
        self.assertIsNotNone(lo.non_minimal_solver)

    def test_od_checks_nm_solver_not_primary(self):
        with self.assertRaises(ValueError):
            homlib.LOMSAC(
                homlib.Wadenback3DV2026OneSided(),
                non_minimal_solver=homlib.Wadenback3DV2026OneSided(),
                non_minimal_settings={"strategy": "overdetermined"},
                inlier_threshold=1.0)

    def test_nm_solver_with_settings_no_strategy_auto_sets_od(self):
        lo = homlib.LOMSAC(
            homlib.Wadenback3DV2026OneSided(),
            non_minimal_solver=homlib.FitzgibbonCVPR2001OneSided(),
            non_minimal_settings={"refine": False},
            inlier_threshold=1.0)
        self.assertEqual(lo.non_minimal_strategy,
                         homlib.NonMinimalStrategy.Overdetermined)

    def test_nm_solver_with_cpp_settings_default_strategy_raises(self):
        nm = homlib.NonMinimalSettings()
        nm.refine = False
        with self.assertRaises(ValueError, msg="NMS object strategy=Sample"):
            homlib.LOMSAC(
                homlib.Wadenback3DV2026OneSided(),
                non_minimal_solver=homlib.FitzgibbonCVPR2001OneSided(),
                non_minimal_settings=nm,
                inlier_threshold=1.0)

    def test_nm_solver_with_cpp_settings_explicit_od_ok(self):
        nm = homlib.NonMinimalSettings()
        nm.strategy = homlib.NonMinimalStrategy.Overdetermined
        nm.refine = False
        lo = homlib.LOMSAC(
            homlib.Wadenback3DV2026OneSided(),
            non_minimal_solver=homlib.FitzgibbonCVPR2001OneSided(),
            non_minimal_settings=nm,
            inlier_threshold=1.0)
        self.assertEqual(lo.non_minimal_strategy,
                         homlib.NonMinimalStrategy.Overdetermined)

    def test_nm_solver_with_explicit_od_and_refine_off(self):
        lo = homlib.LOMSAC(
            homlib.Wadenback3DV2026OneSided(),
            non_minimal_solver=homlib.FitzgibbonCVPR2001OneSided(),
            non_minimal_settings={"strategy": "overdetermined", "refine": False},
            inlier_threshold=1.0)
        self.assertEqual(lo.non_minimal_strategy,
                         homlib.NonMinimalStrategy.Overdetermined)
        self.assertFalse(lo._nm_settings.refine)


@unittest.skipUnless(homlib._HAS_GCRANSAC, "GC-RANSAC not built")
class GCRANSACValidationTestCase(unittest.TestCase):
    """Same validation rules apply to GC-RANSAC."""

    def test_od_with_non_od_solver_raises(self):
        with self.assertRaises(ValueError):
            homlib.GCRANSAC(homlib.Wadenback3DV2026OneSided(),
                            non_minimal_settings={"strategy": "overdetermined"},
                            threshold=1.0)

    def test_nm_solver_with_sample_raises(self):
        with self.assertRaises(ValueError):
            homlib.GCRANSAC(
                homlib.Wadenback3DV2026OneSided(),
                non_minimal_solver=homlib.FitzgibbonCVPR2001OneSided(),
                non_minimal_settings={"strategy": "sample"},
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
                          non_minimal_settings={"strategy": "overdetermined"},
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
                              non_minimal_settings={"strategy": "overdetermined"},
                              inlier_threshold=1.0)

    def test_nm_solver_auto_sets_od(self):
        sa = homlib.SupeRANSAC(
            homlib.Wadenback3DV2026OneSided(),
            non_minimal_solver=homlib.FitzgibbonCVPR2001OneSided(),
            inlier_threshold=1.0)
        self.assertEqual(sa.non_minimal_strategy,
                         homlib.NonMinimalStrategy.Overdetermined)


class NmSolverWiringTestCase(unittest.TestCase):
    """Verify nm_solver is actually wired through to affine/orientation backends."""

    def _make_problem(self, distortion_case=None, n_points=80, n_outliers=15):
        if distortion_case is None:
            distortion_case = homlib.DistortionCase.ONE_SIDED_RIGHT
        config = homlib.ProblemConfig(distortion_case, 0.0, n_points)
        inst = homlib.generate_problem_instance(config)
        x = np.array(inst.x1).T
        y = np.array(inst.x2).T
        A = np.array([np.asarray(a).flatten() for a in inst.A]).T
        ori = np.column_stack(inst.ori)
        rng = np.random.RandomState(42)
        x[:, :n_outliers] = rng.randn(2, n_outliers) * 5
        y[:, :n_outliers] = rng.randn(2, n_outliers) * 5
        return x, y, A, ori, inst

    def test_lomsac_affine_barath_with_nm_solver(self):
        x, y, A, _, inst = self._make_problem(
            distortion_case=homlib.DistortionCase.NO_DISTORTION)
        solver = homlib.BarathVISAPP2016Affine()
        nm = homlib.DLT()
        lo = homlib.LOMSAC(solver, non_minimal_solver=nm,
                           inlier_threshold=1e-3)
        est, stats = lo.estimate(x, y, A=A)
        self.assertGreater(len(stats.inlier_indices), 30)

    def test_lomsac_affine_vo2026_with_nm_solver(self):
        x, y, A, _, inst = self._make_problem()
        solver = homlib.ValtonenOrnhagICPR2026AffineOneSidedRight()
        nm = homlib.FitzgibbonCVPR2001OneSidedRight()
        lo = homlib.LOMSAC(solver, non_minimal_solver=nm,
                           inlier_threshold=1e-3)
        est, stats = lo.estimate(x, y, A=A)
        self.assertGreater(len(stats.inlier_indices), 30)

    def test_lomsac_orientation_vo2026_with_nm_solver(self):
        x, y, _, ori, inst = self._make_problem()
        solver = homlib.ValtonenOrnhagICPR2026OrientationOneSidedRight()
        nm = homlib.FitzgibbonCVPR2001OneSidedRight()
        lo = homlib.LOMSAC(solver, non_minimal_solver=nm,
                           inlier_threshold=1e-3)
        est, stats = lo.estimate(x, y, ori=ori)
        self.assertGreater(len(stats.inlier_indices), 30)

    @unittest.skipUnless(homlib._HAS_GCRANSAC, "GC-RANSAC not built")
    def test_gcransac_affine_with_nm_solver(self):
        x, y, A, _, inst = self._make_problem(
            distortion_case=homlib.DistortionCase.NO_DISTORTION)
        solver = homlib.BarathVISAPP2016Affine()
        nm = homlib.DLT()
        gc = homlib.GCRANSAC(solver, non_minimal_solver=nm, threshold=1e-3)
        est, stats = gc.estimate(x, y, A=A)
        self.assertGreater(len(stats.inliers), 30)

    @unittest.skipUnless(homlib._HAS_GCRANSAC, "GC-RANSAC not built")
    def test_gcransac_orientation_with_nm_solver(self):
        x, y, _, ori, inst = self._make_problem()
        solver = homlib.ValtonenOrnhagICPR2026OrientationOneSidedRight()
        nm = homlib.FitzgibbonCVPR2001OneSidedRight()
        gc = homlib.GCRANSAC(solver, non_minimal_solver=nm, threshold=1e-3)
        est, stats = gc.estimate(x, y, ori=ori)
        self.assertGreater(len(stats.inliers), 30)

    @unittest.skipUnless(homlib._HAS_SUPERANSAC, "SupeRANSAC not built")
    def test_superansac_affine_with_nm_solver(self):
        x, y, A, _, inst = self._make_problem(
            distortion_case=homlib.DistortionCase.NO_DISTORTION)
        solver = homlib.BarathVISAPP2016Affine()
        nm = homlib.DLT()
        sa = homlib.SupeRANSAC(solver, non_minimal_solver=nm,
                               inlier_threshold=1e-3)
        est, stats = sa.estimate(x, y, A=A)
        self.assertGreater(len(stats.inliers), 30)

    @unittest.skipUnless(homlib._HAS_MAGSAC, "MAGSAC not built")
    def test_magsac_affine_with_nm_solver(self):
        x, y, A, _, inst = self._make_problem(
            distortion_case=homlib.DistortionCase.NO_DISTORTION)
        solver = homlib.BarathVISAPP2016Affine()
        nm = homlib.DLT()
        mg = homlib.MAGSAC(solver, non_minimal_solver=nm, sigma_max=1e-3)
        est, stats = mg.estimate(x, y, A=A)
        self.assertGreater(len(stats.inliers), 30)

    # --- Same-type nm_solver tests (affine+affine, ori+ori) ---

    def test_lomsac_barath_affine_with_barath_affine_nm(self):
        """Barath affine as both minimal and non-minimal solver."""
        x, y, A, _, inst = self._make_problem(
            distortion_case=homlib.DistortionCase.NO_DISTORTION)
        solver = homlib.BarathVISAPP2016Affine()
        nm = homlib.BarathVISAPP2016Affine()
        lo = homlib.LOMSAC(solver, non_minimal_solver=nm,
                           inlier_threshold=1e-3)
        est, stats = lo.estimate(x, y, A=A)
        self.assertGreater(len(stats.inlier_indices), 30)

    def test_lomsac_vo_affine_with_vo_affine_nm(self):
        """VÖ affine as both minimal and non-minimal solver."""
        x, y, A, _, inst = self._make_problem()
        solver = homlib.ValtonenOrnhagICPR2026AffineOneSidedRight()
        nm = homlib.ValtonenOrnhagICPR2026AffineOneSidedRight()
        lo = homlib.LOMSAC(solver, non_minimal_solver=nm,
                           inlier_threshold=1e-3)
        est, stats = lo.estimate(x, y, A=A)
        self.assertGreater(len(stats.inlier_indices), 30)

    def test_lomsac_vo_orientation_with_vo_orientation_nm(self):
        """VÖ orientation as both minimal and non-minimal solver."""
        x, y, _, ori, inst = self._make_problem()
        solver = homlib.ValtonenOrnhagICPR2026OrientationOneSidedRight()
        nm = homlib.ValtonenOrnhagICPR2026OrientationOneSidedRight()
        lo = homlib.LOMSAC(solver, non_minimal_solver=nm,
                           inlier_threshold=1e-3)
        est, stats = lo.estimate(x, y, ori=ori)
        self.assertGreater(len(stats.inlier_indices), 30)

    def test_lomsac_vo_affine_with_nakano_nm(self):
        """VÖ affine + Nakano (point-based, matching ONE_SIDED_RIGHT)."""
        x, y, A, _, inst = self._make_problem()
        solver = homlib.ValtonenOrnhagICPR2026AffineOneSidedRight()
        nm = homlib.NakanoICPR2025OneSidedRight()
        lo = homlib.LOMSAC(solver, non_minimal_solver=nm,
                           inlier_threshold=1e-3)
        est, stats = lo.estimate(x, y, A=A)
        self.assertGreater(len(stats.inlier_indices), 30)

    def test_lomsac_vo_orientation_with_nakano_nm(self):
        """VÖ orientation + Nakano (point-based, matching ONE_SIDED_RIGHT)."""
        x, y, _, ori, inst = self._make_problem()
        solver = homlib.ValtonenOrnhagICPR2026OrientationOneSidedRight()
        nm = homlib.NakanoICPR2025OneSidedRight()
        lo = homlib.LOMSAC(solver, non_minimal_solver=nm,
                           inlier_threshold=1e-3)
        est, stats = lo.estimate(x, y, ori=ori)
        self.assertGreater(len(stats.inlier_indices), 30)

    def test_vo_with_dlt_nm_distortion_mismatch_raises(self):
        """VÖ (ONE_SIDED_RIGHT) + DLT (NO_DISTORTION) must be rejected."""
        with self.assertRaises(ValueError, msg="distortion_case mismatch"):
            homlib.LOMSAC(
                homlib.ValtonenOrnhagICPR2026AffineOneSidedRight(),
                non_minimal_solver=homlib.DLT(),
                inlier_threshold=1.0)

    def test_vo_with_wadenback_nm_no_od_raises(self):
        """Wadenbäck doesn't support overdetermined — rejected as nm_solver."""
        with self.assertRaises(ValueError, msg="no OD support"):
            homlib.LOMSAC(
                homlib.ValtonenOrnhagICPR2026AffineOneSidedRight(),
                non_minimal_solver=homlib.Wadenback3DV2026OneSidedRight(),
                inlier_threshold=1.0)


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


class KRangePropagationTestCase(unittest.TestCase):
    """k_range propagation from main solver to non-minimal solver."""

    def test_string_nm_solver_inherits_k_range(self):
        s = homlib.FitzgibbonCVPR2001OneSided()
        s.set_k_range(-0.5, 0.05)
        lo = homlib.LOMSAC(s, non_minimal_solver="FitzgibbonCVPR2001OneSided",
                           inlier_threshold=1.0)
        self.assertEqual(lo.non_minimal_solver.k_range, (-0.5, 0.05))

    def test_solver_id_nm_solver_inherits_k_range(self):
        s = homlib.FitzgibbonCVPR2001OneSided()
        s.set_k_range(-0.3, 0.08)
        lo = homlib.LOMSAC(
            s, non_minimal_solver=homlib.SolverId.FitzgibbonCVPR2001OneSided,
            inlier_threshold=1.0)
        self.assertEqual(lo.non_minimal_solver.k_range, (-0.3, 0.08))

    def test_preconstructed_nm_solver_keeps_own_k_range(self):
        s = homlib.FitzgibbonCVPR2001OneSided()
        s.set_k_range(-0.5, 0.05)
        nm = homlib.FitzgibbonCVPR2001OneSided()
        nm.set_k_range(-0.2, 0.02)
        lo = homlib.LOMSAC(s, non_minimal_solver=nm, inlier_threshold=1.0)
        self.assertEqual(lo.non_minimal_solver.k_range, (-0.2, 0.02))

    def test_no_k_range_no_propagation(self):
        s = homlib.FitzgibbonCVPR2001OneSided()
        lo = homlib.LOMSAC(s, non_minimal_solver="FitzgibbonCVPR2001OneSided",
                           inlier_threshold=1.0)
        self.assertIsNone(lo.non_minimal_solver.k_range)
