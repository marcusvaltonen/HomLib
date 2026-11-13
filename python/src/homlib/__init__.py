from __future__ import annotations

import warnings
from abc import ABC, abstractmethod

from . import _core
from ._core import (
    __doc__,
    __version__,
    PoseData,
    LORansacOptions,
    RansacStatistics,
    DistortionCase,
    NonMinimalStrategy,
    NonMinimalSettings,
    ProblemConfig,
    ProblemInstance,
    generate_problem_instance,
    SolverId,
    FeatureType,
    Solver as _CppSolver,
    dispatch_estimate_point as _dispatch_estimate_point,
    dispatch_estimate_affine as _dispatch_estimate_affine,
    dispatch_estimate_orientation as _dispatch_estimate_orientation,
    dispatch_lomsac_point as _dispatch_lomsac_point,
    dispatch_lomsac_affine as _dispatch_lomsac_affine,
    dispatch_lomsac_orientation as _dispatch_lomsac_orientation,
)

# pybind11 enums share a metaclass that lacks __iter__; patch it once so
# `for val in SomeEnum` works on any pybind11 enum in this package.
_pb11_meta = type(_core.SolverId)
if not hasattr(_pb11_meta, '__iter__'):
    _pb11_meta.__iter__ = lambda cls: iter(cls.__members__.values())

# ---------------------------------------------------------------------------
# Solver base class (Python-side, adds dispatch-based estimate())
# ---------------------------------------------------------------------------


class Solver(_CppSolver):
    """Base class for all HomLib solvers.

    Dispatches :meth:`estimate` based on :attr:`feature_type`
    (``FeatureType.Points``, ``FeatureType.Affine``, or
    ``FeatureType.Orientation``).
    """

    def __init__(self, solver_id, **kwargs):
        if isinstance(solver_id, str):
            member = SolverId.__members__.get(solver_id)
            if member is None:
                raise ValueError(
                    f"Unknown solver {solver_id!r}; "
                    f"expected one of {sorted(SolverId.__members__)}")
            solver_id = member
        super().__init__(solver_id, **kwargs)

    def estimate(self, x, y, A=None, ori=None):
        ft = self.feature_type
        if ft == FeatureType.Points:
            return _dispatch_estimate_point(self, x, y)
        elif ft == FeatureType.Affine:
            if A is None:
                raise ValueError("Affine solver requires A parameter")
            return _dispatch_estimate_affine(self, x, y, A)
        elif ft == FeatureType.Orientation:
            if ori is None:
                raise ValueError("Orientation solver requires ori parameter")
            return _dispatch_estimate_orientation(self, x, y, ori)
        else:
            raise ValueError(f"Unknown feature type: {ft}")


def solver_ids(feature_type=None, names=None):
    """The solvers of HomLib, as a list of :class:`SolverId`, in the order of the enum.

    Parameters
    ----------
    feature_type : FeatureType, optional
        Only the solvers of this feature type, e.g. ``FeatureType.Points`` for the solvers
        that take point matches alone. Default: every solver.
    names : list of str, optional
        Only the solvers with these names (e.g. ``["DLT", "GuoPRL2019"]``, from a command line).
        Raises ``ValueError`` for a name that isn't a solver.

    Examples
    --------
    >>> for solver_id in homlib.solver_ids(homlib.FeatureType.Points):
    ...     lomsac = homlib.LOMSAC(homlib.Solver(solver_id), inlier_threshold=3.0)
    """
    if names is not None:
        unknown = sorted(set(names) - set(SolverId.__members__))
        if unknown:
            raise ValueError(f"unknown solvers: {', '.join(unknown)}")
    return [s for name, s in SolverId.__members__.items()
            if (feature_type is None or Solver(s).feature_type == feature_type)
            and (names is None or name in names)]


# ---------------------------------------------------------------------------
# Point-based solvers
# ---------------------------------------------------------------------------

class CaiPAMI2025(Solver):
    """4-point homography solver using the ACA decomposition by Cai et al. (PAMI 2025). No distortion."""

    def __init__(self):
        super().__init__(SolverId.CaiPAMI2025)


class DLT(Solver):
    """Standard 4-point DLT homography solver. No distortion.

    Parameters
    ----------
    use_qr : bool
        If True, compute the null space via QR decomposition of A^T
        (faster). If False (default), use SVD (more robust).
    """

    def __init__(self, use_qr=False):
        super().__init__(SolverId.DLT, use_qr=use_qr)


class FitzgibbonCVPR2001OneSided(Solver):
    """One-sided radial distortion solver by Fitzgibbon (CVPR 2001), 5+ points."""

    def __init__(self):
        super().__init__(SolverId.FitzgibbonCVPR2001OneSided)


class FitzgibbonCVPR2001TwoSidedEqual(Solver):
    """Two-sided equal radial distortion solver by Fitzgibbon (CVPR 2001), 5+ points."""

    def __init__(self):
        super().__init__(SolverId.FitzgibbonCVPR2001TwoSidedEqual)


class KukelovaCVPR2015TwoSidedEqual(Solver):
    """Two-sided equal radial distortion solver by Kukelova et al. (CVPR 2015), 5-point."""

    def __init__(self):
        super().__init__(SolverId.KukelovaCVPR2015TwoSidedEqual)


class KukelovaCVPR2015TwoSidedEqual6Pt(Solver):
    """Two-sided equal radial distortion solver by Kukelova et al. (CVPR 2015), 6-point."""

    def __init__(self):
        super().__init__(SolverId.KukelovaCVPR2015TwoSidedEqual6Pt)


class KukelovaCVPR2015TwoSided(Solver):
    """Two-sided radial distortion solver by Kukelova et al. (CVPR 2015), 5-point."""

    def __init__(self):
        super().__init__(SolverId.KukelovaCVPR2015TwoSided)


class KukelovaCVPR2015TwoSided6Pt(Solver):
    """Two-sided radial distortion solver by Kukelova et al. (CVPR 2015), 6-point."""

    def __init__(self):
        super().__init__(SolverId.KukelovaCVPR2015TwoSided6Pt)


class NakanoICPR2025OneSided(Solver):
    """One-sided radial distortion solver by Nakano (ICPR 2025), 5+ points."""

    def __init__(self, extra_check=False):
        super().__init__(SolverId.NakanoICPR2025OneSided, extra_check=extra_check)


class NakanoICPR2025OneSidedRight(Solver):
    """One-sided right radial distortion solver by Nakano (ICPR 2025), 5+ points.

    Automatically derived from :class:`NakanoICPR2025OneSided` by swapping the
    point arguments and inverting the homography.
    """

    def __init__(self, extra_check=False):
        super().__init__(SolverId.NakanoICPR2025OneSidedRight, extra_check=extra_check)


class FitzgibbonCVPR2001OneSidedRight(Solver):
    """One-sided right radial distortion solver by Fitzgibbon (CVPR 2001), 5+ points.

    Automatically derived from :class:`FitzgibbonCVPR2001OneSided` by swapping the
    point arguments and inverting the homography.
    """

    def __init__(self):
        super().__init__(SolverId.FitzgibbonCVPR2001OneSidedRight)


class Wadenback3DV2026OneSided(Solver):
    """One-sided radial distortion solver by Wadenback et al. (3DV 2026), 5-point."""

    def __init__(self, extra_check=True):
        super().__init__(SolverId.Wadenback3DV2026OneSided, extra_check=extra_check)


class Wadenback3DV2026OneSidedRight(Solver):
    """One-sided right radial distortion solver by Wadenback et al. (3DV 2026), 5-point.

    Automatically derived from :class:`Wadenback3DV2026OneSided` by swapping the
    point arguments and inverting the homography.
    """

    def __init__(self, extra_check=True):
        super().__init__(SolverId.Wadenback3DV2026OneSidedRight, extra_check=extra_check)


class Wadenback3DV2026TwoSidedEqual(Solver):
    """Two-sided equal radial distortion solver by Wadenback et al. (3DV 2026), 5-point."""

    def __init__(self, extra_check=True):
        super().__init__(SolverId.Wadenback3DV2026TwoSidedEqual, extra_check=extra_check)


class Wadenback3DV2026TwoSided(Solver):
    """Two-sided radial distortion solver by Wadenback et al. (3DV 2026), 5-point."""

    def __init__(self, extra_check=False):
        super().__init__(SolverId.Wadenback3DV2026TwoSided, extra_check=extra_check)


class GuoPRL2019(Solver):
    """4-point homography solver using Guo's explicit method (PRL 2019). No distortion."""

    def __init__(self):
        super().__init__(SolverId.GuoPRL2019)


# ---------------------------------------------------------------------------
# Affine-covariant solvers
# ---------------------------------------------------------------------------

class BarathVISAPP2016Affine(Solver):
    """Affine homography solver by Barath and Hajder (VISAPP 2016), 2-point."""

    def __init__(self):
        super().__init__(SolverId.BarathVISAPP2016Affine)


class ValtonenOrnhagICPR2026AffineOneSided(Solver):
    """Affine-covariant one-sided radial distortion solver by Valtonen Ornhag (ICPR 2026), 2-point.

    Distortion is on the left (y) side. Automatically derived from the native
    right-sided solver by swapping the point arguments, inverting the affine
    features and homography.
    """

    def __init__(self, extra_check=False):
        super().__init__(SolverId.ValtonenOrnhagICPR2026AffineOneSided, extra_check=extra_check)


class ValtonenOrnhagICPR2026AffineOneSidedRight(Solver):
    """Affine-covariant one-sided right radial distortion solver by Valtonen Ornhag (ICPR 2026), 2-point.

    Distortion is on the right (x) side. This is the native solver direction.
    """

    def __init__(self, extra_check=False):
        super().__init__(SolverId.ValtonenOrnhagICPR2026AffineOneSidedRight, extra_check=extra_check)


# ---------------------------------------------------------------------------
# Orientation-covariant solvers
# ---------------------------------------------------------------------------

class ValtonenOrnhagICPR2026OrientationOneSided(Solver):
    """Orientation-covariant one-sided radial distortion solver by Valtonen Ornhag (ICPR 2026), 4-point.

    Distortion is on the left (y) side. Automatically derived from the native
    right-sided solver by swapping the point arguments, the orientation
    components, and inverting the homography.
    """

    def __init__(self):
        super().__init__(SolverId.ValtonenOrnhagICPR2026OrientationOneSided)


class ValtonenOrnhagICPR2026OrientationOneSidedRight(Solver):
    """Orientation-covariant one-sided right radial distortion solver by Valtonen Ornhag (ICPR 2026), 4-point.

    Distortion is on the right (x) side. This is the native solver direction.
    """

    def __init__(self):
        super().__init__(SolverId.ValtonenOrnhagICPR2026OrientationOneSidedRight)


# ---------------------------------------------------------------------------
# RANSAC backends
# ---------------------------------------------------------------------------

class RobustEstimator(ABC):
    """Abstract base class for RANSAC-style robust estimation backends.

    All backends share a common interface: wrap a solver at construction,
    then call :meth:`estimate` with point correspondences (and optional
    affine/orientation data).

    Subclasses must implement :meth:`estimate`.
    """

    def __init__(self, solver, *, non_minimal_solver=None,
                 non_minimal_settings=None):
        self._solver = solver if isinstance(solver, Solver) else Solver(solver)
        self._non_minimal_solver, self._nm_settings = _init_non_minimal(
            self._solver, non_minimal_solver, non_minimal_settings)

    @property
    def solver(self):
        """The underlying solver."""
        return self._solver

    @property
    def non_minimal_solver(self):
        """The non-minimal solver, or None."""
        return self._non_minimal_solver

    @property
    def non_minimal_strategy(self):
        """The non-minimal solving strategy."""
        return self._nm_settings.strategy

    @non_minimal_strategy.setter
    def non_minimal_strategy(self, value):
        parsed = _parse_strategy(value)
        old = self._nm_settings.strategy
        self._nm_settings.strategy = parsed
        try:
            _check_non_minimal_params(
                self._solver, self._non_minimal_solver, self._nm_settings)
        except ValueError:
            self._nm_settings.strategy = old
            raise

    @property
    def non_minimal_settings(self):
        """The non-minimal settings."""
        return self._nm_settings

    @abstractmethod
    def estimate(self, x, y, A=None, ori=None, scores=None):
        """Run robust estimation.

        Parameters
        ----------
        x : numpy.ndarray
            2xN array of points in the first image.
        y : numpy.ndarray
            2xN array of points in the second image.
        A : numpy.ndarray, optional
            4xN local affine transformations (required for affine solvers).
        ori : numpy.ndarray, optional
            2xN orientation vectors (required for orientation solvers).
        scores : numpy.ndarray, optional
            1-D array of per-match confidence scores in [0, 1],
            sorted in descending order.  Used by score-aware
            samplers (PROSAC, ARSampler, ImportanceSampler).
            Supported by GC-RANSAC, MAGSAC, and SupeRANSAC.

        Returns
        -------
        tuple[PoseData, object]
            Best model and backend-specific estimation statistics.
        """


_NON_MINIMAL_STRATEGY_STRINGS = {
    "sample": NonMinimalStrategy.Sample,
    "overdetermined": NonMinimalStrategy.Overdetermined,
}


def _parse_strategy(value):
    """Parse a strategy value (string or enum) to NonMinimalStrategy."""
    if isinstance(value, NonMinimalStrategy):
        return value
    if isinstance(value, str):
        key = value.lower()
        if key not in _NON_MINIMAL_STRATEGY_STRINGS:
            raise ValueError(
                f"Unknown strategy {value!r}; "
                f"expected one of {list(_NON_MINIMAL_STRATEGY_STRINGS)}"
            )
        return _NON_MINIMAL_STRATEGY_STRINGS[key]
    raise TypeError(
        f"strategy must be a str or NonMinimalStrategy, "
        f"got {type(value).__name__}"
    )


def _resolve_non_minimal_settings(value):
    """Resolve a NonMinimalSettings, dict, or None to a NonMinimalSettings."""
    if value is None:
        return NonMinimalSettings()
    if isinstance(value, NonMinimalSettings):
        return NonMinimalSettings(strategy=value.strategy,
                                  sample_size=value.sample_size,
                                  refine=value.refine)
    if isinstance(value, dict):
        kw = dict(value)
        if "strategy" in kw:
            kw["strategy"] = _parse_strategy(kw["strategy"])
        return NonMinimalSettings(**kw)
    raise TypeError(
        f"non_minimal_settings must be a NonMinimalSettings or dict, "
        f"got {type(value).__name__}"
    )


def _resolve_non_minimal_solver(value, minimal_solver):
    """Resolve a non-minimal solver specification to a Solver instance."""
    if value is None:
        return None
    if isinstance(value, (str, SolverId)):
        value = Solver(value)
        if minimal_solver.k_range is not None:
            value.set_k_range(*minimal_solver.k_range)
    if not isinstance(value, Solver):
        raise TypeError(
            f"non_minimal_solver must be a Solver, SolverId, or str, "
            f"got {type(value).__name__}"
        )
    if (value.feature_type != FeatureType.Points
            and value.feature_type != minimal_solver.feature_type):
        raise ValueError(
            f"non_minimal_solver feature_type ({value.feature_type.name}) "
            f"must be Points or match the minimal solver's "
            f"({minimal_solver.feature_type.name})"
        )
    if value.distortion_case != minimal_solver.distortion_case:
        raise ValueError(
            f"non_minimal_solver distortion_case "
            f"({value.distortion_case.name}) must match the minimal "
            f"solver's ({minimal_solver.distortion_case.name})"
        )
    if not value.supports_overdetermined:
        raise ValueError(
            f"{value.id.name} does not support overdetermined solving; "
            f"choose a solver with supports_overdetermined=True"
        )
    return value


def _check_non_minimal_params(solver, non_minimal_solver, nm_settings):
    """Validate the combination of non-minimal solver and strategy."""
    strategy = nm_settings.strategy
    if (non_minimal_solver is not None
            and strategy != NonMinimalStrategy.Overdetermined):
        raise ValueError(
            "non_minimal_solver requires strategy='overdetermined'"
        )
    if (strategy == NonMinimalStrategy.Overdetermined
            and non_minimal_solver is None
            and not solver.supports_overdetermined):
        raise ValueError(
            f"{solver.id.name} does not support overdetermined solving; "
            f"use strategy='sample' or a solver that accepts N > minimal"
        )


def _init_non_minimal(solver, non_minimal_solver, non_minimal_settings):
    """Resolve, auto-set strategy, and validate non-minimal parameters."""
    nm_solver = _resolve_non_minimal_solver(non_minimal_solver, solver)
    nm_settings = _resolve_non_minimal_settings(non_minimal_settings)
    strategy_explicit = (
        isinstance(non_minimal_settings, NonMinimalSettings)
        or (isinstance(non_minimal_settings, dict)
            and "strategy" in non_minimal_settings))
    if nm_solver is not None and not strategy_explicit:
        nm_settings.strategy = NonMinimalStrategy.Overdetermined
    _check_non_minimal_params(solver, nm_solver, nm_settings)
    return nm_solver, nm_settings


def _warn_sample_size_ignored(nm_settings, backend_name):
    if nm_settings.sample_size > 0:
        warnings.warn(
            f"NonMinimalSettings.sample_size is ignored by {backend_name}; "
            f"only LOMSAC honours this setting",
            stacklevel=3,
        )


class LOMSAC(RobustEstimator):
    """Locally Optimized MSAC (LOMSAC) robust estimation backend.

    Parameters
    ----------
    solver : Solver or str or SolverId
        A solver instance, name, or enum (any feature type).
    non_minimal_solver : Solver or str or SolverId, optional
        A separate solver for non-minimal (overdetermined) estimation.
        Must be a point solver or match the primary solver's feature type
        (affine or orientation), and have the same ``distortion_case``.
        Sets strategy to ``'overdetermined'`` when the settings don't
        specify one.
    non_minimal_settings : NonMinimalSettings or dict, optional
        Settings for the non-minimal step: ``strategy`` (``"sample"`` or
        ``"overdetermined"``), ``sample_size``, ``refine``.
    inlier_threshold : float, optional
        Convenience parameter: raw inlier threshold (not squared).
        Internally squared and stored as ``squared_inlier_threshold``.
        Mutually exclusive with ``squared_inlier_threshold``.
    options : LORansacOptions, optional
        Pre-configured options object. If not given, a default is created.
    **kwargs
        Any attribute of :class:`LORansacOptions` can be passed as a keyword
        argument to override defaults (e.g. ``squared_inlier_threshold=0.01``).
    """

    NonMinimalStrategy = NonMinimalStrategy

    def __init__(self, solver, *, non_minimal_solver=None,
                 non_minimal_settings=None, options=None, **kwargs):
        super().__init__(solver, non_minimal_solver=non_minimal_solver,
                         non_minimal_settings=non_minimal_settings)
        self._options = options if options is not None else LORansacOptions()
        if "inlier_threshold" in kwargs:
            if "squared_inlier_threshold" in kwargs:
                raise ValueError(
                    "Cannot specify both inlier_threshold and "
                    "squared_inlier_threshold"
                )
            t = kwargs.pop("inlier_threshold")
            self._options.squared_inlier_threshold = t * t
        for key, value in kwargs.items():
            if not hasattr(self._options, key):
                raise ValueError(f"Unknown LOMSAC option: {key}")
            setattr(self._options, key, value)

    @property
    def options(self):
        """The LOMSAC options."""
        return self._options

    def estimate(self, x, y, A=None, ori=None, scores=None):
        ft = self._solver.feature_type
        if ft == FeatureType.Points:
            return _dispatch_lomsac_point(
                self._solver, x, y, self._nm_settings,
                self._non_minimal_solver, self._options)
        elif ft == FeatureType.Affine:
            if A is None:
                raise ValueError("Affine solver requires A parameter")
            return _dispatch_lomsac_affine(
                self._solver, x, y, A, self._nm_settings,
                self._non_minimal_solver, self._options)
        elif ft == FeatureType.Orientation:
            if ori is None:
                raise ValueError("Orientation solver requires ori parameter")
            return _dispatch_lomsac_orientation(
                self._solver, x, y, ori, self._nm_settings,
                self._non_minimal_solver, self._options)
        else:
            raise ValueError(f"Unknown feature type: {ft}")

    def __repr__(self):
        return f"LOMSAC(solver={self._solver!r})"


# ---------------------------------------------------------------------------
# SupeRANSAC backend (optional)
# ---------------------------------------------------------------------------

_HAS_SUPERANSAC = getattr(_core, "_has_superansac", False)

if _HAS_SUPERANSAC:
    from ._core import (
        SuperansacSettings,
        SuperansacStatistics,
        SuperansacScoringType,
        SuperansacSamplerType,
        SuperansacLocalOptType,
        SuperansacNeighborhoodType,
        SuperansacNeighborhoodSettings,
        SuperansacTerminationType,
        SuperansacInlierSelectorType,
        SuperansacARSamplerSettings,
        SuperansacLOSettings,
        SuperansacDegensacSettings,
        dispatch_superansac_point as _dispatch_superansac_point,
        dispatch_superansac_affine as _dispatch_superansac_affine,
        dispatch_superansac_orientation as _dispatch_superansac_orientation,
    )


class SupeRANSAC(RobustEstimator):
    """SupeRANSAC robust estimation backend.

    Requires HomLib to be built with ``-DHOMLIB_USE_SUPERANSAC=ON``.

    Parameters
    ----------
    solver : Solver or str or SolverId
        A solver instance, name, or enum (any feature type).
    non_minimal_solver : Solver or str or SolverId, optional
        A separate solver for non-minimal (overdetermined) estimation.
        Must be a point solver or match the primary solver's feature type
        (affine or orientation), and have the same ``distortion_case``.
        Sets strategy to ``'overdetermined'`` when the settings don't
        specify one.
    non_minimal_settings : NonMinimalSettings or dict, optional
        Settings for the non-minimal step: ``strategy`` (``"sample"`` or
        ``"overdetermined"``), ``sample_size``, ``refine``.
        Note: ``sample_size`` is ignored by SupeRANSAC (only LOMSAC honours it).
    settings : SuperansacSettings, optional
        Pre-configured settings object. If not given, a default is created.
    **kwargs
        Any attribute of :class:`SuperansacSettings` can be passed as a keyword
        argument to override defaults (e.g. ``inlier_threshold=1.5``).
    """

    if _HAS_SUPERANSAC:
        Settings = SuperansacSettings
        Statistics = SuperansacStatistics
        ScoringType = SuperansacScoringType
        SamplerType = SuperansacSamplerType
        LocalOptType = SuperansacLocalOptType
        NeighborhoodType = SuperansacNeighborhoodType
        NeighborhoodSettings = SuperansacNeighborhoodSettings
        TerminationType = SuperansacTerminationType
        InlierSelectorType = SuperansacInlierSelectorType
        ARSamplerSettings = SuperansacARSamplerSettings
        LOSettings = SuperansacLOSettings
        DegensacSettings = SuperansacDegensacSettings

    def __init__(self, solver, *, non_minimal_solver=None,
                 non_minimal_settings=None, settings=None, **kwargs):
        if not _HAS_SUPERANSAC:
            raise ImportError(
                "HomLib was built without SupeRANSAC support. "
                "Rebuild with -DHOMLIB_USE_SUPERANSAC=ON"
            )
        super().__init__(solver, non_minimal_solver=non_minimal_solver,
                         non_minimal_settings=non_minimal_settings)
        _warn_sample_size_ignored(self._nm_settings, "SupeRANSAC")
        self._settings = settings if settings is not None else SuperansacSettings()
        for key, value in kwargs.items():
            if not hasattr(self._settings, key):
                raise ValueError(f"Unknown SupeRANSAC setting: {key}")
            setattr(self._settings, key, value)

    @property
    def settings(self):
        """The SupeRANSAC settings."""
        return self._settings

    def estimate(self, x, y, A=None, ori=None, scores=None):
        probs = list(scores) if scores is not None else []
        ft = self._solver.feature_type
        if ft == FeatureType.Points:
            return _dispatch_superansac_point(
                self._solver, x, y, probs, self._nm_settings,
                self._non_minimal_solver, self._settings)
        elif ft == FeatureType.Affine:
            if A is None:
                raise ValueError("Affine solver requires A parameter")
            return _dispatch_superansac_affine(
                self._solver, x, y, A, probs, self._nm_settings,
                self._non_minimal_solver, self._settings)
        elif ft == FeatureType.Orientation:
            if ori is None:
                raise ValueError("Orientation solver requires ori parameter")
            return _dispatch_superansac_orientation(
                self._solver, x, y, ori, probs, self._nm_settings,
                self._non_minimal_solver, self._settings)
        else:
            raise ValueError(f"Unknown feature type: {ft}")

    def __repr__(self):
        return f"SupeRANSAC(solver={self._solver!r})"


# ---------------------------------------------------------------------------
# GC-RANSAC backend (optional)
# ---------------------------------------------------------------------------

_HAS_GCRANSAC = getattr(_core, "_has_gcransac", False)

if _HAS_GCRANSAC:
    from ._core import (
        GCRANSACSettings,
        GCRANSACStatistics,
        GCSamplerType,
        dispatch_gcransac_point as _dispatch_gcransac_point,
        dispatch_gcransac_affine as _dispatch_gcransac_affine,
        dispatch_gcransac_orientation as _dispatch_gcransac_orientation,
    )


class GCRANSAC(RobustEstimator):
    """GC-RANSAC robust estimation backend.

    Requires HomLib to be built with ``-DHOMLIB_USE_GCRANSAC=ON``.

    Parameters
    ----------
    solver : Solver or str or SolverId
        A solver instance, name, or enum (any feature type).
    non_minimal_solver : Solver or str or SolverId, optional
        A separate solver for non-minimal (overdetermined) estimation.
        Must be a point solver or match the primary solver's feature type
        (affine or orientation), and have the same ``distortion_case``.
        Sets strategy to ``'overdetermined'`` when the settings don't
        specify one.
    non_minimal_settings : NonMinimalSettings or dict, optional
        Settings for the non-minimal step: ``strategy`` (``"sample"`` or
        ``"overdetermined"``), ``sample_size``, ``refine``.
        Note: ``sample_size`` is ignored by GC-RANSAC (only LOMSAC honours it).
    settings : GCRANSACSettings, optional
        Pre-configured settings object. If not given, a default is created.
    **kwargs
        Any attribute of :class:`GCRANSACSettings` can be passed as a keyword
        argument to override defaults (e.g. ``threshold=3.0``, or ``seed=0`` for
        reproducible results).
    """

    if _HAS_GCRANSAC:
        Settings = GCRANSACSettings
        Statistics = GCRANSACStatistics
        SamplerType = GCSamplerType
        NonMinimalStrategy = NonMinimalStrategy

    def __init__(self, solver, *, non_minimal_solver=None,
                 non_minimal_settings=None, settings=None, **kwargs):
        if not _HAS_GCRANSAC:
            raise ImportError(
                "HomLib was built without GC-RANSAC support. "
                "Rebuild with -DHOMLIB_USE_GCRANSAC=ON"
            )
        super().__init__(solver, non_minimal_solver=non_minimal_solver,
                         non_minimal_settings=non_minimal_settings)
        _warn_sample_size_ignored(self._nm_settings, "GC-RANSAC")
        self._settings = settings if settings is not None else GCRANSACSettings()
        for key, value in kwargs.items():
            if not hasattr(self._settings, key):
                raise ValueError(f"Unknown GC-RANSAC setting: {key}")
            setattr(self._settings, key, value)

    @property
    def settings(self):
        """The GC-RANSAC settings."""
        return self._settings

    def estimate(self, x, y, A=None, ori=None, scores=None):
        probs = list(scores) if scores is not None else []
        ft = self._solver.feature_type
        if ft == FeatureType.Points:
            return _dispatch_gcransac_point(
                self._solver, x, y, probs, self._nm_settings,
                self._non_minimal_solver, self._settings)
        elif ft == FeatureType.Affine:
            if A is None:
                raise ValueError("Affine solver requires A parameter")
            return _dispatch_gcransac_affine(
                self._solver, x, y, A, probs, self._nm_settings,
                self._non_minimal_solver, self._settings)
        elif ft == FeatureType.Orientation:
            if ori is None:
                raise ValueError("Orientation solver requires ori parameter")
            return _dispatch_gcransac_orientation(
                self._solver, x, y, ori, probs, self._nm_settings,
                self._non_minimal_solver, self._settings)
        else:
            raise ValueError(f"Unknown feature type: {ft}")

    def __repr__(self):
        return f"GCRANSAC(solver={self._solver!r})"


# ---------------------------------------------------------------------------
# MAGSAC backend (optional)
# ---------------------------------------------------------------------------

_HAS_MAGSAC = getattr(_core, "_has_magsac", False)

if _HAS_MAGSAC:
    from ._core import (
        MAGSACSettings,
        MAGSACStatistics,
        dispatch_magsac_point as _dispatch_magsac_point,
        dispatch_magsac_affine as _dispatch_magsac_affine,
        dispatch_magsac_orientation as _dispatch_magsac_orientation,
    )


class MAGSAC(RobustEstimator):
    """MAGSAC/MAGSAC++ robust estimation backend.

    Requires HomLib to be built with ``-DHOMLIB_USE_MAGSAC=ON``.

    Parameters
    ----------
    solver : Solver or str or SolverId
        A solver instance, name, or enum (any feature type).
    non_minimal_solver : Solver or str or SolverId, optional
        A separate solver for non-minimal (overdetermined) estimation.
        Must be a point solver or match the primary solver's feature type
        (affine or orientation), and have the same ``distortion_case``.
        Sets strategy to ``'overdetermined'`` when the settings don't
        specify one.
    non_minimal_settings : NonMinimalSettings or dict, optional
        Settings for the non-minimal step: ``strategy`` (``"sample"`` or
        ``"overdetermined"``), ``sample_size``, ``refine``.
        Note: ``sample_size`` is ignored by MAGSAC (only LOMSAC honours it).
    settings : MAGSACSettings, optional
        Pre-configured settings object. If not given, a default is created.
    **kwargs
        Any attribute of :class:`MAGSACSettings` can be passed as a keyword
        argument to override defaults (e.g. ``sigma_max=10.0``, or ``seed=0`` for
        reproducible results).
    """

    if _HAS_MAGSAC:
        Settings = MAGSACSettings
        Statistics = MAGSACStatistics
    if _HAS_GCRANSAC:
        SamplerType = GCSamplerType
        NonMinimalStrategy = NonMinimalStrategy

    def __init__(self, solver, *, non_minimal_solver=None,
                 non_minimal_settings=None, settings=None, **kwargs):
        if not _HAS_MAGSAC:
            raise ImportError(
                "HomLib was built without MAGSAC support. "
                "Rebuild with -DHOMLIB_USE_MAGSAC=ON"
            )
        super().__init__(solver, non_minimal_solver=non_minimal_solver,
                         non_minimal_settings=non_minimal_settings)
        _warn_sample_size_ignored(self._nm_settings, "MAGSAC")
        self._settings = settings if settings is not None else MAGSACSettings()
        for key, value in kwargs.items():
            if not hasattr(self._settings, key):
                raise ValueError(f"Unknown MAGSAC setting: {key}")
            setattr(self._settings, key, value)

    @property
    def settings(self):
        """The MAGSAC settings."""
        return self._settings

    def estimate(self, x, y, A=None, ori=None, scores=None):
        probs = list(scores) if scores is not None else []
        ft = self._solver.feature_type
        if ft == FeatureType.Points:
            return _dispatch_magsac_point(
                self._solver, x, y, probs, self._nm_settings,
                self._non_minimal_solver, self._settings)
        elif ft == FeatureType.Affine:
            if A is None:
                raise ValueError("Affine solver requires A parameter")
            return _dispatch_magsac_affine(
                self._solver, x, y, A, probs, self._nm_settings,
                self._non_minimal_solver, self._settings)
        elif ft == FeatureType.Orientation:
            if ori is None:
                raise ValueError("Orientation solver requires ori parameter")
            return _dispatch_magsac_orientation(
                self._solver, x, y, ori, probs, self._nm_settings,
                self._non_minimal_solver, self._settings)
        else:
            raise ValueError(f"Unknown feature type: {ft}")

    def __repr__(self):
        return f"MAGSAC(solver={self._solver!r})"


__all__ = [
    "__doc__",
    "__version__",
    # Data types
    "PoseData",
    "LORansacOptions",
    "RansacStatistics",
    # Enums
    "SolverId",
    "FeatureType",
    "DistortionCase",
    "NonMinimalStrategy",
    "NonMinimalSettings",
    # Base class
    "Solver",
    "solver_ids",
    # Problem generation
    "ProblemConfig",
    "ProblemInstance",
    "generate_problem_instance",
    # Point-based solvers
    "CaiPAMI2025",
    "DLT",
    "FitzgibbonCVPR2001OneSided",
    "FitzgibbonCVPR2001TwoSidedEqual",
    "KukelovaCVPR2015TwoSidedEqual",
    "KukelovaCVPR2015TwoSidedEqual6Pt",
    "KukelovaCVPR2015TwoSided",
    "KukelovaCVPR2015TwoSided6Pt",
    "NakanoICPR2025OneSided",
    "NakanoICPR2025OneSidedRight",
    "FitzgibbonCVPR2001OneSidedRight",
    "Wadenback3DV2026OneSided",
    "Wadenback3DV2026OneSidedRight",
    "Wadenback3DV2026TwoSidedEqual",
    "Wadenback3DV2026TwoSided",
    "GuoPRL2019",
    # Affine-covariant solvers
    "BarathVISAPP2016Affine",
    "ValtonenOrnhagICPR2026AffineOneSided",
    "ValtonenOrnhagICPR2026AffineOneSidedRight",
    # Orientation-covariant solvers
    "ValtonenOrnhagICPR2026OrientationOneSided",
    "ValtonenOrnhagICPR2026OrientationOneSidedRight",
    # RANSAC backends
    "RobustEstimator",
    "LOMSAC",
    "SupeRANSAC",
    "GCRANSAC",
    "MAGSAC",
    # Conditional availability flags
    "_HAS_SUPERANSAC",
    "_HAS_GCRANSAC",
    "_HAS_MAGSAC",
]

if _HAS_SUPERANSAC:
    __all__ += [
        "SuperansacSettings",
        "SuperansacStatistics",
        "SuperansacScoringType",
        "SuperansacSamplerType",
        "SuperansacLocalOptType",
    ]

if _HAS_GCRANSAC:
    __all__ += [
        "GCRANSACSettings",
        "GCRANSACStatistics",
    ]

if _HAS_MAGSAC:
    __all__ += [
        "MAGSACSettings",
        "MAGSACStatistics",
    ]
