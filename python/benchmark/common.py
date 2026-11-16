"""What the benchmarks share: the robust estimators and their samplers, the description of a run
(its command and its environment), the check of the settings of two runs, and the mapping of points
with an estimated model."""

import datetime
import os
import platform
import shlex
import subprocess
import sys
from pathlib import Path

import numpy as np

import homlib


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


# The samplers of each backend, by their names on the command line, and their names in the backend's
# sampler enum (GCSamplerType, SuperansacSamplerType); LOMSAC has only uniform sampling. PROSAC draws
# its first samples from the first matches: it needs them sorted by decreasing score, and doesn't use
# the scores themselves (prepare_matches sorts them). The importance and AR samplers use the scores.
GC_SAMPLERS = {"uniform": "Uniform", "prosac": "PROSAC", "importance": "ImportanceSampler", "ar": "ARSampler"}
SAMPLERS = {"lomsac": {"uniform": None}, "gcransac": GC_SAMPLERS, "magsac": GC_SAMPLERS,
            "superansac": GC_SAMPLERS}
SORTED_SAMPLERS = {"prosac"}   # take the matches sorted by decreasing score
SCORED_SAMPLERS = {"importance", "ar"}   # take the scores


def add_sampler_argument(parser):
    """Adds --sampler to an argparse parser; check_sampler checks it against the backend."""
    parser.add_argument("--sampler", choices=sorted({s for samplers in SAMPLERS.values() for s in samplers}),
                        default="uniform",
                        help="Sampler of the robust estimator (default: uniform, for every backend). prosac gets "
                             "the matches sorted by decreasing score; importance and ar get the scores. LOMSAC "
                             "has only uniform")


def check_sampler(backend, sampler):
    """Raises a ValueError if the backend has no such sampler."""
    if sampler not in SAMPLERS[backend]:
        raise ValueError(f"--backend {backend} has no sampler {sampler} (its samplers: "
                         f"{', '.join(sorted(SAMPLERS[backend]))})")


def prepare_matches(x, y, scores, sampler):
    """The matches (2xN x and y) and the scores to give to estimate() with the sampler, and their
    order: order[i] is the index in x and y of the i-th match returned (None if the order is
    unchanged). PROSAC
    gets the matches sorted by decreasing score (a stable sort: matches of equal scores keep their
    order); only the samplers that use the scores get them (None for the others). original_indices
    maps the inliers back to the indices of x and y. Raises a ValueError if the sampler needs scores
    and there are none."""
    if scores is None and sampler in SORTED_SAMPLERS | SCORED_SAMPLERS:
        raise ValueError(f"the sampler {sampler} needs the scores of the matches, and the matches have none")
    order = None
    if sampler in SORTED_SAMPLERS:
        order = np.argsort(-np.asarray(scores), kind="stable")
        x, y, scores = x[:, order], y[:, order], np.asarray(scores)[order]
    scores = np.ascontiguousarray(scores, dtype=float) if sampler in SCORED_SAMPLERS else None
    return np.ascontiguousarray(x), np.ascontiguousarray(y), scores, order


def with_sampler(settings, superansac="prosac (unsorted)"):
    """The settings of a saved run, with the sampler of the files saved before --sampler: uniform,
    except with SupeRANSAC, whose default sampler is PROSAC: superansac names it, as the benchmark
    gave it the matches (HPatches unsorted, HEB sorted). None if settings is None."""
    if settings is None or "sampler" in settings:
        return settings
    return {**settings, "sampler": superansac if settings.get("backend") == "superansac" else "uniform"}


def original_indices(inliers, order):
    """The indices of the inliers in the matches before prepare_matches reordered them"""
    return np.asarray(inliers, dtype=int) if order is None else order[np.asarray(inliers, dtype=int)]


def make_estimator(backend, solver, threshold, iterations, seed,
                   min_iterations=None, non_minimal_strategy=None,
                   non_minimal_solver=None, sampler="uniform", **extra):
    """The robust estimator, with the inlier threshold, and exactly this many iterations if
    iterations isn't None, or at least min_iterations if it isn't None. sampler: a sampler of
    SAMPLERS[backend] (uniform by default for every backend: SupeRANSAC's own default is PROSAC).
    extra: other settings of the backend, by their names in HomLib."""
    check_sampler(backend, sampler)
    cls, threshold_name, min_name, max_name = BACKENDS[backend]
    options = {threshold_name: threshold, **extra}
    if SAMPLERS[backend][sampler] is not None:
        enum = homlib.SuperansacSamplerType if backend == "superansac" else homlib.GCSamplerType
        options["sampler"] = getattr(enum, SAMPLERS[backend][sampler])
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
