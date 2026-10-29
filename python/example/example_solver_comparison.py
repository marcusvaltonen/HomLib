r"""
Solver comparison
=================

Compare multiple one-sided solvers on synthetic data.
We run 1000 random problem instances and show histograms of
:math:`\log_{10}` homography error and cumulative runtime distributions.

The solvers compared are by Wadenbäck et al. [1]_, Nakano [2]_,
and Fitzgibbon [3]_.

References
----------
.. [1] M. Wadenbäck, M. Valtonen Örnhag, and J. Edstedt, "Radially Distorted
   Homographies, Revisited," in Proc. 3DV, 2026.
.. [2] G. Nakano, "Inverse DLT Method for One-Sided Radial Distortion
   Homography," in Proc. ICPR, 2024.
.. [3] A. Fitzgibbon, "Simultaneous linear estimation of multiple view geometry
   and lens distortion," in Proc. CVPR, 2001.
"""

import time

import numpy as np
import matplotlib.pyplot as plt

import homlib

###############################################################################
# Define solvers to compare

solver_specs = {
    "Wadenbäck (3DV 2026)": lambda x, y: homlib.estimate_wadenback_3dv_2026_one_sided(x, y, True),
    "Nakano (ICPR 2025)": lambda x, y: homlib.estimate_nakano_icpr_2025_one_sided(x, y, True),
    "Fitzgibbon (CVPR 2001)": lambda x, y: homlib.estimate_fitzgibbon_cvpr_2001_one_sided(x, y),
}

###############################################################################
# Run 1000 random instances per solver

n_instances = 1000
all_errors = {name: [] for name in solver_specs}
all_times = {name: [] for name in solver_specs}

for _ in range(n_instances):
    config = homlib.ProblemConfig(homlib.DistortionCase.ONE_SIDED_LEFT, 0.0, 5)
    instance = homlib.generate_problem_instance(config)

    x = np.array(instance.x1).T
    y = np.array(instance.x2).T

    for name, fn in solver_specs.items():
        t0 = time.perf_counter()
        poses = fn(x, y)
        elapsed = (time.perf_counter() - t0) * 1000

        best_err = float("inf")
        for p in poses:
            err = instance.hom_error(p.homography)
            best_err = min(best_err, err)

        all_errors[name].append(best_err)
        all_times[name].append(elapsed)

###############################################################################
# Print summary statistics

print(f"{'Solver':<25s} {'Median H err':>14s} {'Mean time (ms)':>16s}")
print("-" * 58)
for name in solver_specs:
    med_err = np.median(all_errors[name])
    mean_t = np.mean(all_times[name])
    print(f"{name:<25s} {med_err:14.2e} {mean_t:16.4f}")

###############################################################################
# Error distribution
# ------------------
# Histogram of :math:`\log_{10}` homography error across all instances.

colors = ["#4c72b0", "#dd8452", "#55a868"]
bins = np.linspace(-16, 2, 50)

plt.figure(figsize=(7, 4))
for (name, errs), c in zip(all_errors.items(), colors):
    log_errs = np.log10(np.clip(errs, 1e-16, None))
    plt.hist(log_errs, bins=bins, alpha=0.6, label=name, color=c)

plt.xlabel(r"$\log_{10}$(homography error)")
plt.ylabel("Count")
plt.title("Error distribution (1000 instances)")
plt.legend(fontsize=8)
plt.tight_layout()

###############################################################################
# Runtime distribution
# --------------------
# Cumulative distribution of per-instance solver runtime.

plt.figure(figsize=(7, 4))
for (name, times), c in zip(all_times.items(), colors):
    sorted_t = np.sort(times)
    cdf = np.arange(1, len(sorted_t) + 1) / len(sorted_t)
    plt.plot(sorted_t, cdf, label=name, color=c)

plt.xlabel("Runtime (ms)")
plt.ylabel("Cumulative fraction")
plt.title("Runtime distribution (1000 instances)")
plt.legend(fontsize=8)
plt.tight_layout()
