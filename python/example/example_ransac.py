r"""
Synthetic RANSAC experiment
===========================

In this example, we show how to robustly estimate a homography
in the presence of outliers using synthetic data and different
distortion profiles.
"""

import numpy as np

import homlib


solvers = [
    homlib.lomsac_wadenback_3dv_2026_one_sided,
    homlib.lomsac_wadenback_3dv_2026_two_sided_equal,
    homlib.lomsac_wadenback_3dv_2026_two_sided,
]

options = homlib.LORansacOptions()
options.squared_inlier_threshold = 0.005
options.final_least_squares = True

N = 100
H = np.random.randn(3, 3)
x = np.random.randn(2, N)
y = H @ np.vstack((x, np.ones((1, N))))
y = y[:2] / y[2]

print(f"H={H / H[2,2]}")


# Add noise
noise_std = 1e-4
for i in range(N):
    x[:,i] += np.random.randn(2) * noise_std
    y[:,i] += np.random.randn(2) * noise_std


# Add outliers
nbr_outliers = 20
for i in range(nbr_outliers):
    x[:,i] = np.random.randn(2) * 10
    y[:,i] = np.random.randn(2) * 10


for solver in solvers:
    estimate, stats = solver(x, y, options)

    print(f"Results for {solver.__name__}")
    print(estimate)
    print(estimate.homography)
    print(f"H error = {np.linalg.norm(estimate.homography / estimate.homography[2,2] - H / H[2,2]):4e}")
    print(f"Dist. error = {np.linalg.norm(estimate.distortion_parameter):4e}")  # No distortion added - should be zero
    print(stats)
    print(stats.inlier_indices)

