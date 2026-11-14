r"""
Synthetic experiment
====================

In this example, we show how to estimate a homography from synthetic data
using the one-sided solver by Wadenbäck et al. :cite:`Wadenback2026`.
"""

import numpy as np

import homlib

######################################
# Generate synthetic data
N = 5
H = np.random.randn(3, 3)
x = np.random.randn(2, N)
y = H @ np.vstack((x, np.ones((1, N))))
y = y[:2] / y[2]

print(f"H={H / H[2,2]}")

######################################
# Run estimator from homlib on the minimal sample.
solver = homlib.Wadenback3DV2026OneSided()
poses = solver.estimate(x, y)

######################################
# Test against ground truth

for i, p in enumerate(poses):
    print(f"Pose {i}: {p} with H =")
    print(p.homography / p.homography[2, 2])
    print(f"H error = {np.linalg.norm(p.homography / p.homography[2, 2] - H / H[2, 2]):4e}")
    print(f"Dist. error = {np.linalg.norm(p.distortion_parameter2):4e}")  # No distortion added - should be zero
