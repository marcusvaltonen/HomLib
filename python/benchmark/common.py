"""What the benchmarks share: the robust estimators, the check of the settings of two runs, and the
mapping of points with an estimated model."""

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


def make_estimator(backend, solver, threshold, iterations, seed):
    """The robust estimator, with the inlier threshold, and exactly this many iterations if
    iterations isn't None"""
    cls, threshold_name, min_name, max_name = BACKENDS[backend]
    options = {threshold_name: threshold}
    if iterations is not None:
        options.update({min_name: iterations, max_name: iterations})
    if backend == "lomsac":
        options["random_seed"] = seed
    return getattr(homlib, cls)(solver, **options)


def inlier_indices(stats):
    """The indices of the inliers that the estimator found"""
    return stats.inlier_indices if hasattr(stats, "inlier_indices") else stats.inliers


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
