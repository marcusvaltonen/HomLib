r"""
RANSAC backend comparison
=========================

Compare the four robust estimation backends available in HomLib on
synthetic data with known ground truth and one-sided radial distortion.

Thresholds are set to the normalised-coordinate equivalent of 3 pixels
at a simulated focal length of 1000 px, i.e.
:math:`\tau = 5 / f = 0.005`.

The backends are LOMSAC [1]_, SupeRANSAC [2]_, GC-RANSAC [3]_, and
MAGSAC++ [4]_.  LOMSAC is always available; the others require building
with their respective CMake flags (see :ref:`installation`).

References
----------
.. [1] K. Lebeda, J. Matas, and O. Chum, "Fixing the Locally Optimized
   RANSAC," in Proc. BMVC, 2012.
.. [2] D. Barath, "SupeRANSAC: One RANSAC to Rule Them All," arXiv
   preprint arXiv:2506.04803, 2025.
.. [3] D. Barath and J. Matas, "Graph-Cut RANSAC," in Proc. CVPR, 2018.
.. [4] D. Barath, J. Noskova, M. Ivashechkin, and J. Matas, "MAGSAC++, a
   Fast, Reliable and Accurate Robust Estimator," in Proc. CVPR, 2020.
"""

import time
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt

import homlib

np.random.seed(42)

###############################################################################
# Check which optional backends are available.

has_superansac = homlib._HAS_SUPERANSAC
has_gcransac = homlib._HAS_GCRANSAC
has_magsac = homlib._HAS_MAGSAC

print(
    f"Backends: LOMSAC=True, SupeRANSAC={has_superansac}, "
    f"GC-RANSAC={has_gcransac}, MAGSAC++={has_magsac}"
)

###############################################################################
# Configuration
# -------------
# Thresholds are expressed as pixel equivalents at a focal length of
# 1000 px and converted to normalised coordinates.

N_POINTS = 100
N_RUNS = 5
NOISE = 1e-3
OUTLIER_RATIOS = [0.3, 0.5, 0.7]

FOCAL_LENGTH = 1000.0
PIXEL_THRESHOLD = 5.0
THRESHOLD = PIXEL_THRESHOLD / FOCAL_LENGTH

solver = homlib.Wadenback3DV2026OneSided()

###############################################################################
# Create backend instances
# ------------------------
# All backends use the same 3-pixel equivalent threshold.

backends = {
    "LOMSAC": homlib.LOMSAC(
        solver, inlier_threshold=THRESHOLD
    ),
}

if has_superansac:
    backends["SupeRANSAC"] = homlib.SupeRANSAC(
        solver, inlier_threshold=THRESHOLD
    )
if has_gcransac:
    backends["GC-RANSAC"] = homlib.GCRANSAC(
        solver, threshold=THRESHOLD
    )
if has_magsac:
    backends["MAGSAC++"] = homlib.MAGSAC(
        solver, sigma_max=THRESHOLD
    )


###############################################################################
# Helpers
# -------


def generate_data(noise, n_points=N_POINTS):
    """Generate a synthetic problem instance with radial distortion."""
    config = homlib.ProblemConfig(
        homlib.DistortionCase.ONE_SIDED_LEFT, noise, n_points
    )
    instance = homlib.generate_problem_instance(config)
    return {
        "x": np.array(instance.x1).T,
        "y": np.array(instance.x2).T,
        "H_gt": instance.posedata.homography,
        "k_gt": instance.posedata.distortion_parameter2,
    }


def add_outliers(data, ratio):
    """Return a copy of *data* with a fraction *ratio* of outliers."""
    n = data["x"].shape[1]
    n_out = int(n * ratio)
    idx = np.random.choice(n, n_out, replace=False)
    out = {k: (v.copy() if hasattr(v, "copy") else v) for k, v in data.items()}
    out["x"][:, idx] = np.random.randn(2, n_out) * 2
    out["y"][:, idx] = np.random.randn(2, n_out) * 2
    return out


###############################################################################
# Sweep outlier ratios
# --------------------
# Run each backend on the same data and record homography error,
# distortion error, and wall-clock time.

results = defaultdict(lambda: defaultdict(list))

for ratio in OUTLIER_RATIOS:
    data = add_outliers(generate_data(NOISE), ratio)
    H_gt_norm = data["H_gt"] / data["H_gt"][2, 2]

    for name, backend in backends.items():
        errors, dist_errors, times = [], [], []
        for _ in range(N_RUNS):
            t0 = time.perf_counter()
            est, _ = backend.estimate(data["x"], data["y"])
            elapsed = time.perf_counter() - t0

            H_est = est.homography / est.homography[2, 2]
            errors.append(np.linalg.norm(H_est - H_gt_norm))
            dist_errors.append(abs(est.distortion_parameter2 - data["k_gt"]))
            times.append(elapsed * 1000)

        results[name]["ratios"].append(ratio)
        results[name]["errors"].append(np.median(errors))
        results[name]["dist_errors"].append(np.median(dist_errors))
        results[name]["times"].append(np.median(times))

    print(f"  {ratio:.0%} outliers done")

###############################################################################
# Summary table
# -------------

print(
    f"\n{'Backend':<14s} {'Outliers':>8s} {'H err':>10s} "
    f"{'k err':>10s} {'ms':>8s}"
)
print("-" * 54)
for name in backends:
    d = results[name]
    for i, ratio in enumerate(d["ratios"]):
        print(
            f"{name:<14s} {ratio:>7.0%} {d['errors'][i]:>10.2e}"
            f" {d['dist_errors'][i]:>10.2e} {d['times'][i]:>8.1f}"
        )

###############################################################################
# Plots
# -----
# Three panels: homography accuracy, distortion accuracy, and speed as a
# function of outlier ratio.

STYLES = {
    "LOMSAC": {"ls": "-", "marker": "o", "color": "#1f77b4"},
    "SupeRANSAC": {"ls": "--", "marker": "s", "color": "#ff7f0e"},
    "GC-RANSAC": {"ls": "-.", "marker": "^", "color": "#2ca02c"},
    "MAGSAC++": {"ls": ":", "marker": "D", "color": "#d62728"},
}

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

for name, d in results.items():
    s = STYLES[name]
    kw = dict(
        linestyle=s["ls"], marker=s["marker"], color=s["color"],
        label=name, markersize=5,
    )
    axes[0].plot(d["ratios"], d["errors"], **kw)
    axes[1].plot(d["ratios"], d["dist_errors"], **kw)
    axes[2].plot(d["ratios"], d["times"], **kw)

axes[0].set(
    xlabel="Outlier ratio", ylabel="Median H error",
    title="Homography accuracy", yscale="log",
)
axes[1].set(
    xlabel="Outlier ratio", ylabel=r"Median $|k_{est} - k_{gt}|$",
    title="Distortion accuracy", yscale="log",
)
axes[2].set(
    xlabel="Outlier ratio", ylabel="Median time (ms)",
    title="Speed",
)

for ax in axes:
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(fontsize=7)

plt.tight_layout()
plt.show()
