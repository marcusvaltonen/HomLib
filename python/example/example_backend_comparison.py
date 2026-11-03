r"""
RANSAC backend comparison
=========================

Minimal example showing how to use the four robust estimation backends
in HomLib on the same synthetic data.

All backends wrap the same solver and produce a
:class:`~homlib.PoseData` estimate.  LOMSAC is always available; the
others require building with their respective CMake flags.
"""

import numpy as np

import homlib

np.random.seed(42)

###############################################################################
# Generate synthetic data with outliers

solver = homlib.Wadenback3DV2026OneSided()

N = 500
H_gt = np.array([
    [1.2, 0.1, 0.05],
    [0.05, 0.9, -0.02],
    [0.001, -0.001, 1.0],
])

x = np.random.randn(2, N) * 0.3
y_hom = H_gt @ np.vstack((x, np.ones((1, N))))
y = y_hom[:2] / y_hom[2]

x += np.random.randn(2, N) * 1e-4
y += np.random.randn(2, N) * 1e-4

n_outliers = 100
x[:, :n_outliers] = np.random.randn(2, n_outliers) * 5
y[:, :n_outliers] = np.random.randn(2, n_outliers) * 5

H_gt_norm = H_gt / H_gt[2, 2]

###############################################################################
# LOMSAC (always available)
# -------------------------

lomsac = homlib.LOMSAC(solver, inlier_threshold=0.07)
est, stats = lomsac.estimate(x, y)
H_est = est.homography / est.homography[2, 2]
print(f"LOMSAC    — H error: {np.linalg.norm(H_est - H_gt_norm):.2e}, "
      f"inliers: {len(stats.inlier_indices)}")

###############################################################################
# SupeRANSAC (optional)
# ---------------------

if homlib._HAS_SUPERANSAC:
    sa = homlib.SupeRANSAC(solver, inlier_threshold=0.07)
    # PROSAC (default) expects correspondences sorted by quality;
    # use Uniform for unsorted synthetic data.
    sa.settings.sampler = homlib.SuperansacSamplerType.Uniform
    est, stats = sa.estimate(x, y)
    H_est = est.homography / est.homography[2, 2]
    print(f"SupeRANSAC — H error: {np.linalg.norm(H_est - H_gt_norm):.2e}, "
          f"inliers: {len(stats.inliers)}")

###############################################################################
# GC-RANSAC (optional)
# --------------------

if homlib._HAS_GCRANSAC:
    gc = homlib.GCRANSAC(solver, threshold=0.07)
    est, stats = gc.estimate(x, y)
    H_est = est.homography / est.homography[2, 2]
    print(f"GC-RANSAC  — H error: {np.linalg.norm(H_est - H_gt_norm):.2e}, "
          f"inliers: {len(stats.inliers)}")

###############################################################################
# MAGSAC++ (optional)
# -------------------

if homlib._HAS_MAGSAC:
    mg = homlib.MAGSAC(solver, sigma_max=0.5)
    est, stats = mg.estimate(x, y)
    H_est = est.homography / est.homography[2, 2]
    print(f"MAGSAC++   — H error: {np.linalg.norm(H_est - H_gt_norm):.2e}, "
          f"inliers: {len(stats.inliers)}")

###############################################################################
# .. note::
#
#    On clean synthetic data the backends perform similarly, since outliers
#    are easy to separate.  The more sophisticated backends — GC-RANSAC,
#    MAGSAC++, and SupeRANSAC — show their real advantage on real-world
#    images where noise is spatially correlated and the inlier/outlier
#    boundary is less clear-cut.  See the ``example_real_data.py`` gallery
#    example for a comparison on real images.
