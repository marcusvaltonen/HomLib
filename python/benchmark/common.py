"""What the benchmarks share: the robust estimators, and the check of the settings of two runs."""

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


def make_estimator(backend, solver, threshold, iterations, seed, non_minimal_strategy=None):
    """The robust estimator, with the inlier threshold, and exactly this many iterations if
    iterations isn't None"""
    cls, threshold_name, min_name, max_name = BACKENDS[backend]
    options = {threshold_name: threshold}
    if iterations is not None:
        options.update({min_name: iterations, max_name: iterations})
    if backend == "lomsac":
        options["random_seed"] = seed
    if non_minimal_strategy is not None:
        options["non_minimal_strategy"] = non_minimal_strategy
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
