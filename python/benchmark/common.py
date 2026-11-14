"""What the benchmarks share: the robust estimators, the description of a run (its command and its
environment), the check of the settings of two runs, the mapping of points with an estimated
model, the ABCs for feature matchers, and the match result types."""

import datetime
import os
import platform
import shlex
import subprocess
import sys
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path

import numpy as np

import homlib


# ---------------------------------------------------------------------------
# Match result types
# ---------------------------------------------------------------------------

@dataclass
class PointMatchResult:
    """Matched point correspondences between two images."""
    points1: np.ndarray
    points2: np.ndarray
    scores: np.ndarray | None = None


@dataclass(kw_only=True)
class OrientationMatchResult(PointMatchResult):
    """Point matches with per-keypoint orientation (radians)."""
    angles1: np.ndarray
    angles2: np.ndarray


@dataclass(kw_only=True)
class ScaleOrientationMatchResult(OrientationMatchResult):
    """Point matches with per-keypoint orientation and scale."""
    scales1: np.ndarray
    scales2: np.ndarray


@dataclass(kw_only=True)
class AffineMatchResult(ScaleOrientationMatchResult):
    """Point matches with per-keypoint affine frames (2x2)."""
    affines1: np.ndarray
    affines2: np.ndarray

    def relative_affines(self):
        """affines2 @ inv(affines1) per correspondence, as Nx2x2."""
        return self.affines2 @ np.linalg.inv(self.affines1)

    @classmethod
    def from_lafs(cls, laf1, laf2, scores=None):
        """Construct from LAFs (Nx2x3): centroid = col 2, affine = cols 0:2,
        angle = atan2(a21, a11), scale = sqrt(|det(affine)|)."""
        points1 = laf1[:, :, 2]
        points2 = laf2[:, :, 2]
        aff1 = laf1[:, :, :2]
        aff2 = laf2[:, :, :2]
        angles1 = np.arctan2(aff1[:, 1, 0], aff1[:, 0, 0])
        angles2 = np.arctan2(aff2[:, 1, 0], aff2[:, 0, 0])
        scales1 = np.sqrt(np.abs(np.linalg.det(aff1)))
        scales2 = np.sqrt(np.abs(np.linalg.det(aff2)))
        return cls(points1=points1, points2=points2, scores=scores,
                   angles1=angles1, angles2=angles2,
                   scales1=scales1, scales2=scales2,
                   affines1=aff1, affines2=aff2)


# ---------------------------------------------------------------------------
# Feature matcher ABC
# ---------------------------------------------------------------------------

class FeatureMatcher(ABC):
    """A feature type with extraction and matching.  Subclasses are named by feature type
    (SIFT, SuperPoint, AffNet, ...), not by matcher.

    Subclasses that should be available from the CLI declare ``cli_name`` (the ``--features``
    value), ``cli_matcher`` (the ``--matcher`` value, or None for end-to-end matchers that
    ignore it), and a ``from_args(args)`` classmethod that builds an instance from parsed
    CLI arguments.  ``__init_subclass__`` auto-registers them in ``_cli_registry``."""
    supported_features: set
    _cli_registry: dict = {}

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        cli_name = getattr(cls, "cli_name", None)
        if cli_name is not None:
            cli_matcher = getattr(cls, "cli_matcher", None)
            FeatureMatcher._cli_registry[(cli_name, cli_matcher)] = cls

    def __init__(self):
        self._cache = {}

    @abstractmethod
    def _extract(self, path):
        """Extract features from a single image.  Returns a dict of keypoints, descriptors,
        and any extra data the matcher needs."""
        ...

    def extract(self, path):
        """Cached _extract(): holds up to 2 entries and evicts the most-recently inserted
        one when full, so that the first image (the anchor in a 1-vs-many sequence) stays
        cached while the second image rotates."""
        if path not in self._cache:
            if len(self._cache) >= 2:
                self._cache.pop(next(reversed(self._cache)))
            self._cache[path] = self._extract(path)
        return self._cache[path]

    @abstractmethod
    def _match_features(self, features1, features2):
        """Match two sets of extracted features.  Returns a PointMatchResult or subclass.
        Typically delegates to a private matcher method (_match_ratio_test, _match_lightglue,
        etc.)."""
        ...

    def match(self, path1, path2):
        """Extract from both images, then match.  End-to-end systems (LoMa, ALoMa) override
        this entirely."""
        f1 = self.extract(path1)
        f2 = self.extract(path2)
        return self._match_features(f1, f2)

    def __call__(self, path1, path2):
        return self.match(path1, path2)


# ---------------------------------------------------------------------------
# Feature upgrader ABC
# ---------------------------------------------------------------------------

class FeatureUpgrader(ABC):
    """Enriches match results with additional geometric features."""
    input_type: type
    output_type: type

    @abstractmethod
    def upgrade(self, result, **kwargs):
        ...


# The robust estimators: their names in HomLib, and the names of their inlier threshold and of their
# iteration bounds. MAGSAC++ has no threshold, but a largest noise level, sigma_max, which also
# selects its final inliers (the matches within sigma_max). Only LOMSAC takes a random seed: the
# others draw theirs at random, so their runs are not reproducible.
BACKENDS = {
    "lomsac": ("LOMSAC", "inlier_threshold", "min_num_iterations", "max_num_iterations"),
    "gcransac": ("GCRANSAC", "threshold", "min_iterations", "max_iterations"),
    "magsac": ("MAGSAC", "sigma_max", "min_iterations", "max_iterations"),
    "superansac": ("SupeRANSAC", "inlier_threshold", "min_iterations", "max_iterations"),
}


def make_estimator(backend, solver, threshold, iterations, seed,
                   min_iterations=None, non_minimal_strategy=None,
                   non_minimal_solver=None, **extra):
    """The robust estimator, with the inlier threshold, and exactly this many iterations if
    iterations isn't None, or at least min_iterations if it isn't None. extra: other settings of the
    backend, by their names in HomLib."""
    cls, threshold_name, min_name, max_name = BACKENDS[backend]
    options = {threshold_name: threshold, **extra}
    if iterations is not None:
        options.update({min_name: iterations, max_name: iterations})
    elif min_iterations is not None:
        options[min_name] = min_iterations
    if backend == "lomsac":
        options["random_seed"] = seed
    nm_settings = {}
    if non_minimal_strategy is not None:
        nm_settings["strategy"] = non_minimal_strategy
    if nm_settings:
        options["non_minimal_settings"] = nm_settings
    if non_minimal_solver is not None:
        options["non_minimal_solver"] = non_minimal_solver
    return getattr(homlib, cls)(solver, **options)


def inlier_indices(stats):
    """The indices of the inliers that the estimator found"""
    return stats.inlier_indices if hasattr(stats, "inlier_indices") else stats.inliers


def cpu_model():
    """The model of the CPU, or platform.processor() (often only the architecture) if unknown"""
    try:
        if sys.platform == "darwin":
            return subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True,
                                  text=True, check=True).stdout.strip()
        if sys.platform.startswith("linux"):
            with open("/proc/cpuinfo") as f:
                for line in f:
                    if line.startswith("model name"):
                        return line.split(":", 1)[1].strip()
    except (OSError, subprocess.SubprocessError):
        pass
    return platform.processor()


def git_commit():
    """The commit of the repository of the scripts, with "-dirty" if it has uncommitted changes, or
    None outside a git checkout. It can differ from the HomLib that is installed."""
    def git(*args):
        return subprocess.run(["git", "-C", str(Path(__file__).parent), *args], capture_output=True,
                              text=True, check=True).stdout.strip()
    try:
        return git("rev-parse", "HEAD") + ("-dirty" if git("status", "--porcelain", "--untracked-files=no") else "")
    except (OSError, subprocess.SubprocessError):
        return None


def opencv_version():
    try:
        import cv2
        return cv2.__version__
    except ImportError:
        return None


def run_info():
    """The description of a run, saved with its results: the command line, and the environment (the
    host, the versions of HomLib, of its backends and of the Python packages, and the commit of the
    scripts). The environment doesn't count in the check of the settings of two runs."""
    environment = dict(
        date=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        os=platform.platform(), cpu=cpu_model(), cpu_count=os.cpu_count(),
        python=platform.python_version(), homlib=homlib.__version__,
        backends=sorted(name for name, (cls, *_) in BACKENDS.items() if backend_available(cls)),
        numpy=np.__version__, opencv=opencv_version(),
        scripts_commit=git_commit())
    return dict(command=shlex.join([Path(sys.executable).name, *sys.argv]), environment=environment)


def backend_available(cls):
    """Whether HomLib was built with this robust estimator"""
    return cls == "LOMSAC" or bool(getattr(homlib, f"_HAS_{cls.upper()}", False))


def print_runs(run_a, run_b):
    """Prints the commands of two runs, and the differences between their environments (files saved
    before the description of the runs have none)."""
    for label, run in (("A", run_a), ("B", run_b)):
        if run:
            print(f"{label}: {run['command']}")
    if run_a and run_b:
        ea, eb = run_a["environment"], run_b["environment"]
        differ = [f"{k} ({ea.get(k)} vs {eb.get(k)})" for k in sorted(set(ea) | set(eb))
                  if k != "date" and ea.get(k) != eb.get(k)]
        if differ:
            print("The runs have different environments: " + ", ".join(differ))


def warn_settings(settings_a, settings_b):
    """Prints a warning if two runs have different settings, or if one of them has none."""
    if settings_a is None or settings_b is None:
        print("Warning: a file has no settings (saved by an older version of the script): can't check "
              "that the runs have the same settings")
    elif settings_a != settings_b:
        differ = [f"{k} ({settings_a.get(k)} vs {settings_b.get(k)})"
                  for k in sorted(set(settings_a) | set(settings_b)) if settings_a.get(k) != settings_b.get(k)]
        print("Warning: the runs have different settings: " + ", ".join(differ))


def map_points(pose, p):
    """The points p of image 1 mapped to image 2 with an estimated model: undistorted with k1,
    mapped with H, distorted with k2 (division model); and whether each point has an image. It
    hasn't if H maps it to infinity, or if k2 has no distorted point for it (1 - 4 k2 |u|^2 < 0).
    Both are tested before the division and the square root, so that no NaN or infinity is
    created."""
    k1, k2 = pose.distortion_parameter, pose.distortion_parameter2
    # Undistorted with k1 in homogeneous coordinates (w = 1 + k1 |p|^2 as third coordinate): no
    # division by w
    q = np.array(pose.homography) @ np.vstack([p, 1 + k1 * np.sum(p**2, 0)])
    ok = q[2] != 0
    u = np.divide(q[:2], q[2], out=np.zeros_like(q[:2]), where=ok)
    ok &= np.all(np.abs(u) < 1e150, 0)   # a tiny q[2]: |u|^2 would overflow
    u[:, ~ok] = 0
    disc = 1 - 4 * k2 * np.sum(u**2, 0)
    ok &= disc >= 0
    return u * 2 / (1 + np.sqrt(np.where(ok, disc, 0))), ok
