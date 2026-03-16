import numpy as np

import homlib


H = np.eye(3)
p = homlib.PoseData(H, 1.0, -0.001, -0.002)
print(f"PoseData: {p}")

N = 5
H = np.random.randn(3, 3)
x = np.random.randn(2, N)
y = H @ np.vstack((x, np.ones((1, N))))
y = y[:2] / y[2]

print(f"H={H / H[2,2]}")

poses = homlib.estimate_wadenback_2025_one_sided(x, y, True)

for i, p in enumerate(poses):
    print(f"Pose {i}: {p} with H =")
    print(p.homography / p.homography[2,2])
    print(f"H error = {np.linalg.norm(p.homography / p.homography[2,2] - H / H[2,2]):4e}")
    print(f"Dist. error = {np.linalg.norm(p.distortion_parameter):4e}")  # No distortion added - should be zero
