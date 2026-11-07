r"""
.. _sampler_scores_example:

Confidence scores and sampler selection
=======================================

Demonstrates how to pass per-match confidence scores to GC-RANSAC, MAGSAC,
and SupeRANSAC, and how to select the sampling strategy.

When correspondences come from a learned matcher (e.g. LoMa, SuperGlue,
LightGlue), they typically include a confidence score per match.  Sorting
the matches by descending confidence and passing the scores to a
PROSAC-style sampler lets RANSAC focus on the most reliable matches first,
which can improve accuracy — especially when the inlier ratio is low.

Two common score normalization strategies:

- **Real scores**: min-max normalize the raw matcher confidences to [0, 1].
- **Linear scores**: replace the values with ``np.linspace(1.0, 0.01, N)``,
  preserving the sorted order but discarding the magnitudes.  This is the
  approach typically used by Barath et al. in SupeRANSAC / GC-RANSAC.
"""

import numpy as np

import homlib

np.random.seed(42)

###############################################################################
# Generate synthetic data with inliers, outliers, and confidence scores

N_INLIERS = 200
N_OUTLIERS = 80
N = N_INLIERS + N_OUTLIERS

H_gt = np.array([
    [1.2, 0.1, 0.05],
    [0.05, 0.9, -0.02],
    [0.001, -0.001, 1.0],
])

rng = np.random.RandomState(42)
x = rng.randn(2, N) * 0.3
y_hom = H_gt @ np.vstack((x, np.ones((1, N))))
y = y_hom[:2] / y_hom[2]
x += rng.randn(2, N) * 1e-4
y += rng.randn(2, N) * 1e-4

x[:, :N_OUTLIERS] = rng.randn(2, N_OUTLIERS) * 10
y[:, :N_OUTLIERS] = rng.randn(2, N_OUTLIERS) * 10

inlier_mask = np.zeros(N, dtype=bool)
inlier_mask[N_OUTLIERS:] = True
scores_raw = np.where(
    inlier_mask,
    rng.uniform(0.7, 1.0, N),
    rng.uniform(0.0, 0.3, N),
)

order = np.argsort(-scores_raw)
x = np.ascontiguousarray(x[:, order])
y = np.ascontiguousarray(y[:, order])
scores_raw = scores_raw[order]

###############################################################################
# Score normalization
# ~~~~~~~~~~~~~~~~~~~
# Real scores: min-max to [0, 1].
# Linear scores: ``np.linspace(1.0, 0.01, N)`` (Barath et al. convention).

s_min, s_max = scores_raw.min(), scores_raw.max()
scores_real = (scores_raw - s_min) / (s_max - s_min)
scores_linear = np.linspace(1.0, 0.01, N)

solver = homlib.Wadenback3DV2026OneSided()
threshold = 0.07


def h_error(H_est):
    return np.linalg.norm(
        H_est / H_est[2, 2] - H_gt / H_gt[2, 2]
    )


###############################################################################
# GC-RANSAC with different samplers
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if homlib._HAS_GCRANSAC:
    configs = [
        ("Uniform (no scores)", homlib.GCSamplerType.Uniform, None),
        ("PROSAC + real", homlib.GCSamplerType.PROSAC, scores_real),
        ("PROSAC + linear", homlib.GCSamplerType.PROSAC, scores_linear),
        ("ARSampler + real", homlib.GCSamplerType.ARSampler, scores_real),
        ("ImportanceSampler",
         homlib.GCSamplerType.ImportanceSampler, scores_real),
    ]
    for label, sampler, scores in configs:
        gc = homlib.GCRANSAC(
            solver, threshold=threshold, sampler=sampler
        )
        est, stats = gc.estimate(x, y, scores=scores)
        print(
            f"GC-RANSAC  {label:<25s}  "
            f"inliers={len(stats.inliers):>3d}  "
            f"H_err={h_error(est.homography):.4f}"
        )

###############################################################################
# MAGSAC with PROSAC
# ~~~~~~~~~~~~~~~~~~

if homlib._HAS_MAGSAC:
    mg = homlib.MAGSAC(
        solver,
        sigma_max=threshold,
        sampler=homlib.GCSamplerType.PROSAC,
    )
    est, stats = mg.estimate(x, y, scores=scores_real)
    print(
        f"MAGSAC     PROSAC + real               "
        f"inliers={len(stats.inliers):>3d}  "
        f"H_err={h_error(est.homography):.4f}"
    )

###############################################################################
# SupeRANSAC with ARSampler
# ~~~~~~~~~~~~~~~~~~~~~~~~~~

if homlib._HAS_SUPERANSAC:
    sa = homlib.SupeRANSAC(
        solver,
        inlier_threshold=threshold,
        sampler=homlib.SuperansacSamplerType.ARSampler,
    )
    est, stats = sa.estimate(x, y, scores=scores_real)
    print(
        f"SupeRANSAC ARSampler + real             "
        f"inliers={len(stats.inliers):>3d}  "
        f"H_err={h_error(est.homography):.4f}"
    )

###############################################################################
# Error handling
# ~~~~~~~~~~~~~~
# Selecting a score-dependent sampler without providing scores raises
# a ``RuntimeError`` (GC-RANSAC and MAGSAC).

if homlib._HAS_GCRANSAC:
    gc = homlib.GCRANSAC(
        solver,
        threshold=threshold,
        sampler=homlib.GCSamplerType.ARSampler,
    )
    try:
        gc.estimate(x, y)  # no scores → RuntimeError
    except RuntimeError as e:
        print(f"Expected error: {e}")
