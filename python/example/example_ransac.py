r"""
Synthetic RANSAC experiment
===========================

In this example, we show how to robustly estimate a homography
in the presence of outliers using synthetic data and different
distortion profiles.

The one-sided solver by Wadenbäck et al. :cite:`Wadenback2026` is used,
embedded in a LOMSAC framework :cite:`Lebeda2012`.
"""

import numpy as np

import homlib


###############################################################################
# Create a solver and configure LOMSAC

solver = homlib.Wadenback3DV2026OneSided()
lomsac = homlib.LOMSAC(solver, inlier_threshold=0.07, final_least_squares=True)

###############################################################################
# Generate random problem instance (no radial distortion added for simplicity)

N = 100
H = np.random.randn(3, 3)
x = np.random.randn(2, N)
y = H @ np.vstack((x, np.ones((1, N))))
y = y[:2] / y[2]


###############################################################################
# Add noise

noise_std = 1e-4
for i in range(N):
    x[:, i] += np.random.randn(2) * noise_std
    y[:, i] += np.random.randn(2) * noise_std


###############################################################################
# Add outliers

nbr_outliers = 20
for i in range(nbr_outliers):
    x[:, i] = np.random.randn(2) * 10
    y[:, i] = np.random.randn(2) * 10

###############################################################################
# Estimate homography and distortion coefficient

estimate, stats = lomsac.estimate(x, y)

print(f"H error = {np.linalg.norm(estimate.homography / estimate.homography[2, 2] - H / H[2, 2]):4e}")
print(f"Dist. error = {np.linalg.norm(estimate.distortion_parameter2):4e}")  # No distortion added - should be zero
print(stats.inlier_indices)
