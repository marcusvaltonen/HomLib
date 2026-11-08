from __future__ import annotations

from abc import ABC, abstractmethod

from . import _core
from ._core import (
    __doc__,
    __version__,
    PoseData,
    LORansacOptions,
    RansacStatistics,
    DistortionCase,
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

# ---------------------------------------------------------------------------
# Solver base class (Python-side, adds dispatch-based estimate())
# ---------------------------------------------------------------------------


class Solver(_CppSolver):
    """Base class for all HomLib solvers.

    Dispatches :meth:`estimate` based on :attr:`feature_type`
    (``FeatureType.Points``, ``FeatureType.Affine``, or
    ``FeatureType.Orientation``).
    """

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
    """One-sided radial distortion solver by Fitzgibbon (CVPR 2001), 5-point."""

    def __init__(self):
        super().__init__(SolverId.FitzgibbonCVPR2001OneSided)


class FitzgibbonCVPR2001TwoSidedEqual(Solver):
    """Two-sided equal radial distortion solver by Fitzgibbon (CVPR 2001), 5-point."""

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
    """One-sided radial distortion solver by Nakano (ICPR 2025), 5-point."""

    def __init__(self, extra_check=False):
        super().__init__(SolverId.NakanoICPR2025OneSided, extra_check=extra_check)


class NakanoICPR2025OneSidedRight(Solver):
    """One-sided right radial distortion solver by Nakano (ICPR 2025), 5-point."""

    def __init__(self, extra_check=False):
        super().__init__(SolverId.NakanoICPR2025OneSidedRight, extra_check=extra_check)


class Wadenback3DV2026OneSided(Solver):
    """One-sided radial distortion solver by Wadenback et al. (3DV 2026), 5-point."""

    def __init__(self, extra_check=True):
        super().__init__(SolverId.Wadenback3DV2026OneSided, extra_check=extra_check)


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


class ValtonenOrnhagICPR2026AffineSingleSided(Solver):
    """Affine-covariant one-sided radial distortion solver by Valtonen Ornhag (ICPR 2026), 2-point."""

    def __init__(self, extra_check=False):
        super().__init__(SolverId.ValtonenOrnhagICPR2026AffineSingleSided, extra_check=extra_check)


# ---------------------------------------------------------------------------
# Orientation-covariant solvers
# ---------------------------------------------------------------------------

class ValtonenOrnhagICPR2026OrientationSingleSided(Solver):
    """Orientation-covariant one-sided radial distortion solver by Valtonen Ornhag (ICPR 2026), 4-point."""

    def __init__(self):
        super().__init__(SolverId.ValtonenOrnhagICPR2026OrientationSingleSided)


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

    def __init__(self, solver):
        self._solver = solver

    @property
    def solver(self):
        """The underlying solver."""
        return self._solver

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


class LOMSAC(RobustEstimator):
    """Locally Optimized MSAC (LOMSAC) robust estimation backend.

    Parameters
    ----------
    solver : Solver
        A solver instance (any feature type).
    options : LORansacOptions, optional
        Pre-configured options object. If not given, a default is created.
    inlier_threshold : float, optional
        Convenience parameter: raw inlier threshold (not squared).
        Internally squared and stored as ``squared_inlier_threshold``.
        Mutually exclusive with ``squared_inlier_threshold``.
    **kwargs
        Any attribute of :class:`LORansacOptions` can be passed as a keyword
        argument to override defaults (e.g. ``squared_inlier_threshold=0.01``).
    """

    def __init__(self, solver, options=None, **kwargs):
        super().__init__(solver)
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
            return _dispatch_lomsac_point(self._solver, x, y, self._options)
        elif ft == FeatureType.Affine:
            if A is None:
                raise ValueError("Affine solver requires A parameter")
            return _dispatch_lomsac_affine(self._solver, x, y, A, self._options)
        elif ft == FeatureType.Orientation:
            if ori is None:
                raise ValueError("Orientation solver requires ori parameter")
            return _dispatch_lomsac_orientation(self._solver, x, y, ori, self._options)
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
    solver : Solver
        A solver instance (any feature type).
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

    def __init__(self, solver, settings=None, **kwargs):
        if not _HAS_SUPERANSAC:
            raise ImportError(
                "HomLib was built without SupeRANSAC support. "
                "Rebuild with -DHOMLIB_USE_SUPERANSAC=ON"
            )
        super().__init__(solver)
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
            return _dispatch_superansac_point(self._solver, x, y, probs, self._settings)
        elif ft == FeatureType.Affine:
            if A is None:
                raise ValueError("Affine solver requires A parameter")
            return _dispatch_superansac_affine(self._solver, x, y, A, probs, self._settings)
        elif ft == FeatureType.Orientation:
            if ori is None:
                raise ValueError("Orientation solver requires ori parameter")
            return _dispatch_superansac_orientation(self._solver, x, y, ori, probs, self._settings)
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
        NonMinimalStrategy,
        dispatch_gcransac_point as _dispatch_gcransac_point,
        dispatch_gcransac_affine as _dispatch_gcransac_affine,
        dispatch_gcransac_orientation as _dispatch_gcransac_orientation,
    )


class GCRANSAC(RobustEstimator):
    """GC-RANSAC robust estimation backend.

    Requires HomLib to be built with ``-DHOMLIB_USE_GCRANSAC=ON``.

    Parameters
    ----------
    solver : Solver
        A solver instance (any feature type).
    settings : GCRANSACSettings, optional
        Pre-configured settings object. If not given, a default is created.
    **kwargs
        Any attribute of :class:`GCRANSACSettings` can be passed as a keyword
        argument to override defaults (e.g. ``threshold=3.0``).
    """

    if _HAS_GCRANSAC:
        Settings = GCRANSACSettings
        Statistics = GCRANSACStatistics
        SamplerType = GCSamplerType
        NonMinimalStrategy = NonMinimalStrategy

    def __init__(self, solver, settings=None, **kwargs):
        if not _HAS_GCRANSAC:
            raise ImportError(
                "HomLib was built without GC-RANSAC support. "
                "Rebuild with -DHOMLIB_USE_GCRANSAC=ON"
            )
        super().__init__(solver)
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
            return _dispatch_gcransac_point(self._solver, x, y, probs, self._settings)
        elif ft == FeatureType.Affine:
            if A is None:
                raise ValueError("Affine solver requires A parameter")
            return _dispatch_gcransac_affine(self._solver, x, y, A, probs, self._settings)
        elif ft == FeatureType.Orientation:
            if ori is None:
                raise ValueError("Orientation solver requires ori parameter")
            return _dispatch_gcransac_orientation(self._solver, x, y, ori, probs, self._settings)
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
    solver : Solver
        A solver instance (any feature type).
    settings : MAGSACSettings, optional
        Pre-configured settings object. If not given, a default is created.
    **kwargs
        Any attribute of :class:`MAGSACSettings` can be passed as a keyword
        argument to override defaults (e.g. ``sigma_max=10.0``).
    """

    if _HAS_MAGSAC:
        Settings = MAGSACSettings
        Statistics = MAGSACStatistics
    if _HAS_GCRANSAC:
        SamplerType = GCSamplerType
        NonMinimalStrategy = NonMinimalStrategy

    def __init__(self, solver, settings=None, **kwargs):
        if not _HAS_MAGSAC:
            raise ImportError(
                "HomLib was built without MAGSAC support. "
                "Rebuild with -DHOMLIB_USE_MAGSAC=ON"
            )
        super().__init__(solver)
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
            return _dispatch_magsac_point(self._solver, x, y, probs, self._settings)
        elif ft == FeatureType.Affine:
            if A is None:
                raise ValueError("Affine solver requires A parameter")
            return _dispatch_magsac_affine(self._solver, x, y, A, probs, self._settings)
        elif ft == FeatureType.Orientation:
            if ori is None:
                raise ValueError("Orientation solver requires ori parameter")
            return _dispatch_magsac_orientation(self._solver, x, y, ori, probs, self._settings)
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
    "Wadenback3DV2026OneSided",
    "Wadenback3DV2026TwoSidedEqual",
    "Wadenback3DV2026TwoSided",
    "GuoPRL2019",
    # Affine-covariant solvers
    "BarathVISAPP2016Affine",
    "ValtonenOrnhagICPR2026AffineSingleSided",
    # Orientation-covariant solvers
    "ValtonenOrnhagICPR2026OrientationSingleSided",
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
