r"""
LOMSAC vs SupeRANSAC
====================

Compare the two robust estimation backends available in HomLib on synthetic
data with 20% outliers. LOMSAC [1]_ is always available; SupeRANSAC [2]_
requires building with ``-DHOMLIB_USE_SUPERANSAC=ON``.

References
----------
.. [1] K. Lebeda, J. Matas, and O. Chum, "Fixing the Locally Optimized RANSAC,"
   in Proc. BMVC, 2012.
.. [2] D. Barath, "SupeRANSAC: One RANSAC to Rule Them All,"
   arXiv preprint arXiv:2506.04803, 2025.
"""

import time

import numpy as np
import matplotlib.pyplot as plt

import homlib

###############################################################################
# Generate synthetic data
# -----------------------
# A known homography with 100 inliers and 20 outliers.

np.random.seed(42)

N = 100
H_gt = np.array([[1.2, 0.1, 0.05],
                  [0.05, 0.9, -0.02],
                  [0.001, -0.001, 1.0]])

x = np.random.randn(2, N) * 0.3
y_hom = H_gt @ np.vstack((x, np.ones((1, N))))
y = y_hom[:2] / y_hom[2]

noise_std = 1e-4
x += np.random.randn(2, N) * noise_std
y += np.random.randn(2, N) * noise_std

n_outliers = 20
x[:, :n_outliers] = np.random.randn(2, n_outliers) * 10
y[:, :n_outliers] = np.random.randn(2, n_outliers) * 10

H_gt_norm = H_gt / H_gt[2, 2]

###############################################################################
# Define solvers and backends
# ---------------------------
# We test three point-based solvers with both LOMSAC and (optionally)
# SupeRANSAC.

solvers = {
    "Wadenbäck (3DV 2026)": homlib.Wadenback3DV2026OneSided(),
    "Nakano (ICPR 2025)": homlib.NakanoICPR2025OneSided(),
    "Fitzgibbon (CVPR 2001)": homlib.FitzgibbonCVPR2001OneSided(),
}

has_superansac = homlib._HAS_SUPERANSAC

###############################################################################
# Run benchmarks
# --------------
# Each solver+backend pair runs 10 trials. We record homography error,
# number of inliers, and wall-clock time.

n_runs = 10
results = {}

for solver_name, solver in solvers.items():
    lomsac = homlib.LOMSAC(solver, squared_inlier_threshold=0.005)

    lomsac_errors, lomsac_inliers, lomsac_times = [], [], []
    for _ in range(n_runs):
        t0 = time.perf_counter()
        est, stats = lomsac.estimate(x, y)
        elapsed = time.perf_counter() - t0
        H_est = est.homography / est.homography[2, 2]
        lomsac_errors.append(np.linalg.norm(H_est - H_gt_norm))
        lomsac_inliers.append(len(stats.inlier_indices))
        lomsac_times.append(elapsed * 1000)
    results[f"{solver_name} / LOMSAC"] = {
        "errors": lomsac_errors,
        "inliers": lomsac_inliers,
        "times": lomsac_times,
    }

    if has_superansac:
        sa = homlib.SupeRANSAC(
            solver,
            inlier_threshold=0.07,
            sampler=homlib.SuperansacSamplerType.Uniform,
            scoring=homlib.SuperansacScoringType.MSAC,
            local_optimization=homlib.SuperansacLocalOptType.IRLS,
            final_optimization=homlib.SuperansacLocalOptType.IRLS,
        )

        sa_errors, sa_inliers, sa_times = [], [], []
        for _ in range(n_runs):
            t0 = time.perf_counter()
            est, stats = sa.estimate(x, y)
            elapsed = time.perf_counter() - t0
            H_est = est.homography / est.homography[2, 2]
            sa_errors.append(np.linalg.norm(H_est - H_gt_norm))
            sa_inliers.append(len(stats.inliers))
            sa_times.append(elapsed * 1000)
        results[f"{solver_name} / SupeRANSAC"] = {
            "errors": sa_errors,
            "inliers": sa_inliers,
            "times": sa_times,
        }

###############################################################################
# Summary table

print(f"{'Method':<45s} {'Med. H err':>10s} {'Med. inl.':>10s} {'Med. ms':>10s}")
print("-" * 78)
for key, data in results.items():
    print(f"{key:<45s} {np.median(data['errors']):10.2e}"
          f" {np.median(data['inliers']):10.0f}"
          f" {np.median(data['times']):10.2f}")

###############################################################################
# Accuracy comparison
# -------------------
# Median homography error per solver and backend.

solver_names = list(solvers.keys())
lomsac_med_err = [np.median(results[f"{s} / LOMSAC"]["errors"]) for s in solver_names]

x_pos = np.arange(len(solver_names))
width = 0.35

fig, ax = plt.subplots(figsize=(9, 4))
bars_l = ax.bar(x_pos - width / 2, lomsac_med_err, width, label="LOMSAC")

if has_superansac:
    sa_med_err = [np.median(results[f"{s} / SupeRANSAC"]["errors"]) for s in solver_names]
    bars_s = ax.bar(x_pos + width / 2, sa_med_err, width, label="SupeRANSAC")

ax.set_ylabel("Median homography error")
ax.set_title("LOMSAC vs SupeRANSAC — accuracy")
ax.set_xticks(x_pos)
ax.set_xticklabels(solver_names, fontsize=8)
ax.set_yscale("log")
ax.legend()
ax.grid(True, axis="y", alpha=0.3, which="both")
plt.tight_layout()

###############################################################################
# Timing comparison
# -----------------
# Median estimation time per solver and backend.

lomsac_med_time = [np.median(results[f"{s} / LOMSAC"]["times"]) for s in solver_names]

fig, ax = plt.subplots(figsize=(9, 4))
ax.bar(x_pos - width / 2, lomsac_med_time, width, label="LOMSAC")

if has_superansac:
    sa_med_time = [np.median(results[f"{s} / SupeRANSAC"]["times"]) for s in solver_names]
    ax.bar(x_pos + width / 2, sa_med_time, width, label="SupeRANSAC")

ax.set_ylabel("Median time (ms)")
ax.set_title("LOMSAC vs SupeRANSAC — estimation time")
ax.set_xticks(x_pos)
ax.set_xticklabels(solver_names, fontsize=8)
ax.legend()
ax.grid(True, axis="y", alpha=0.3)
plt.tight_layout()
