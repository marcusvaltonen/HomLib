r"""
.. _sampler_scores_example:

Confidence scores and sampler selection
=======================================

Shows how to pass per-match confidence scores to a RANSAC backend
so it samples the most reliable matches first (PROSAC).

We use the Graffiti image pair with SIFT features.  The SNN ratio
from Lowe's ratio test serves as the confidence score.
"""

import cv2
import numpy as np

import homlib

###############################################################################
# Detect SIFT features and match with ratio test

img1 = cv2.imread('img/grafA.png', cv2.IMREAD_GRAYSCALE)
img2 = cv2.imread('img/grafB.png', cv2.IMREAD_GRAYSCALE)

sift = cv2.SIFT_create(8000)
kps1, des1 = sift.detectAndCompute(img1, None)
kps2, des2 = sift.detectAndCompute(img2, None)

matches = cv2.BFMatcher().knnMatch(des1, des2, k=2)

pts1, pts2, scores = [], [], []
for m, n in matches:
    if m.distance < 0.8 * n.distance:
        pts1.append(kps1[m.queryIdx].pt)
        pts2.append(kps2[m.trainIdx].pt)
        scores.append(1.0 - m.distance / n.distance)

x = np.array(pts1).T
y = np.array(pts2).T
scores = np.array(scores)

###############################################################################
# Sort by confidence (descending) and normalize scores to [0, 1]

order = np.argsort(-scores)
x = np.ascontiguousarray(x[:, order])
y = np.ascontiguousarray(y[:, order])
scores = scores[order]

s_min, s_max = scores.min(), scores.max()
scores = (scores - s_min) / (s_max - s_min)

###############################################################################
# Centre the coordinates (homlib expects centred input)

cx, cy = img1.shape[1] / 2, img1.shape[0] / 2
x[0] -= cx
x[1] -= cy
y[0] -= cx
y[1] -= cy

print(f"{x.shape[1]} matches after ratio test")

###############################################################################
# Compare Uniform vs PROSAC sampling with GC-RANSAC

solver = homlib.CaiPAMI2025()

for sampler, label in [
    (homlib.GCSamplerType.Uniform, "Uniform"),
    (homlib.GCSamplerType.PROSAC, "PROSAC"),
]:
    gc = homlib.GCRANSAC(
        solver, threshold=3.0, sampler=sampler)
    est, stats = gc.estimate(x, y, scores=scores)
    print(f"GC-RANSAC {label:<8s}  "
          f"inliers={len(stats.inliers):>4d}  "
          f"iters={stats.num_iterations}")
