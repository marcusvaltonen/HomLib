"""
Example: using confidence scores with different RANSAC samplers.

Demonstrates how to pass per-match confidence scores to GC-RANSAC,
MAGSAC, and SupeRANSAC, and how to select the sampling strategy.

When correspondences come from a learned matcher (e.g. LoMa, SuperGlue,
LightGlue), they typically include a confidence score per match. Sorting
the matches by descending confidence and passing the scores to a
PROSAC-style sampler lets RANSAC focus on the most reliable matches
first, which can improve accuracy -- especially when the inlier ratio
is low.

Two common score normalization strategies:
  - Real scores: min-max normalize the raw matcher confidences to [0,1].
  - Linear scores: replace the values with np.linspace(1.0, 0.01, N),
    preserving the sorted order but discarding the magnitudes. This is
    the approach typically used by Barath et al. in SupeRANSAC/GC-RANSAC.
"""

import numpy as np
import homlib


def make_scene(seed=42, n_inliers=200, n_outliers=80):
    """Generate a synthetic scene with confidence scores."""
    rng = np.random.RandomState(seed)
    H_gt = np.array([
        [1.2, 0.1, 0.05],
        [0.05, 0.9, -0.02],
        [0.001, -0.001, 1.0],
    ])
    N = n_inliers + n_outliers

    x = rng.randn(2, N) * 0.3
    y_hom = H_gt @ np.vstack((x, np.ones((1, N))))
    y = y_hom[:2] / y_hom[2]
    x += rng.randn(2, N) * 1e-4
    y += rng.randn(2, N) * 1e-4

    x[:, :n_outliers] = rng.randn(2, n_outliers) * 10
    y[:, :n_outliers] = rng.randn(2, n_outliers) * 10

    inlier_mask = np.zeros(N, dtype=bool)
    inlier_mask[n_outliers:] = True
    scores = np.where(
        inlier_mask,
        rng.uniform(0.7, 1.0, N),
        rng.uniform(0.0, 0.3, N),
    )

    order = np.argsort(-scores)
    x = np.ascontiguousarray(x[:, order])
    y = np.ascontiguousarray(y[:, order])
    scores = scores[order]

    return x, y, scores, H_gt


def normalize_scores(scores):
    """Min-max normalize to [0, 1]."""
    s_min, s_max = scores.min(), scores.max()
    if s_max > s_min:
        return (scores - s_min) / (s_max - s_min)
    return np.ones_like(scores)


def h_error(H_est, H_gt):
    return np.linalg.norm(
        H_est / H_est[2, 2] - H_gt / H_gt[2, 2])


def main():
    x, y, scores_raw, H_gt = make_scene()
    N = x.shape[1]
    solver = homlib.Wadenback3DV2026OneSided()
    threshold = 0.07

    scores_real = normalize_scores(scores_raw)
    scores_linear = np.linspace(1.0, 0.01, N)

    print(f"Scene: {N} matches, scores in "
          f"[{scores_real.min():.2f}, "
          f"{scores_real.max():.2f}]\n")

    # --- GC-RANSAC with different samplers ---
    if homlib._HAS_GCRANSAC:
        print("=== GC-RANSAC ===")
        GCS = homlib.GCSamplerType
        configs = [
            ("Uniform (no scores)", GCS.Uniform, None),
            ("PROSAC + real", GCS.PROSAC, scores_real),
            ("PROSAC + linear", GCS.PROSAC, scores_linear),
            ("ARSampler + real", GCS.ARSampler, scores_real),
            ("ImportanceSampler",
             GCS.ImportanceSampler, scores_real),
        ]
        for label, sampler, scores in configs:
            gc = homlib.GCRANSAC(
                solver, threshold=threshold,
                sampler=sampler)
            est, stats = gc.estimate(
                x, y, scores=scores)
            err = h_error(est.homography, H_gt)
            n_in = len(stats.inliers)
            print(f"  {label:<25s}  "
                  f"inliers={n_in:>3d}  "
                  f"H_err={err:.4f}")
        print()

    # --- MAGSAC with PROSAC ---
    if homlib._HAS_MAGSAC:
        print("=== MAGSAC ===")
        mg = homlib.MAGSAC(
            solver, sigma_max=threshold,
            sampler=homlib.GCSamplerType.PROSAC)
        est, stats = mg.estimate(
            x, y, scores=scores_real)
        err = h_error(est.homography, H_gt)
        n_in = len(stats.inliers)
        print(f"  PROSAC + real scores  "
              f"inliers={n_in:>3d}  "
              f"H_err={err:.4f}")
        print()

    # --- SupeRANSAC with ARSampler ---
    if homlib._HAS_SUPERANSAC:
        print("=== SupeRANSAC ===")
        sa = homlib.SupeRANSAC(
            solver, inlier_threshold=threshold,
            sampler=homlib.SuperansacSamplerType.ARSampler)
        est, stats = sa.estimate(
            x, y, scores=scores_real)
        err = h_error(est.homography, H_gt)
        n_in = len(stats.inliers)
        print(f"  ARSampler + real      "
              f"inliers={n_in:>3d}  "
              f"H_err={err:.4f}")
        print()

    # --- Error when scores are required but missing ---
    print("=== Error handling ===")
    if homlib._HAS_GCRANSAC:
        gc = homlib.GCRANSAC(
            solver, threshold=threshold,
            sampler=homlib.GCSamplerType.ARSampler)
        try:
            gc.estimate(x, y)
        except RuntimeError as e:
            print(f"  Expected error: {e}")


if __name__ == "__main__":
    main()
