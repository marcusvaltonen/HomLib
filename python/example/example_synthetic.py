r"""
Synthetic experiment
====================

In this example, we show how to estimate a homography from synthetic data.
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
poses = homlib.estimate_wadenback_3dv_2026_one_sided(x, y, True)

######################################
# Test against ground truth

for i, p in enumerate(poses):
    print(f"Pose {i}: {p} with H =")
    print(p.homography / p.homography[2,2])
    print(f"H error = {np.linalg.norm(p.homography / p.homography[2,2] - H / H[2,2]):4e}")
    print(f"Dist. error = {np.linalg.norm(p.distortion_parameter):4e}")  # No distortion added - should be zero
