r"""
Homography from known fundamental matrix (3PT + F)
===================================================

This example demonstrates the **3PT** solver from Barath & Hajder
:cite:`Barath2016` that estimates a homography when the fundamental matrix is
already known.  Given three point correspondences and F, the epipolar
constraint reduces the degrees of freedom so that only three points (instead
of four) are needed.

The same paper also proposes a **HAF** solver that uses a single affine
correspondence together with F — see :class:`~homlib.BarathVISAPP2016HAF`.

.. warning::

   The fundamental matrix is **undefined for purely planar scenes**: when all
   points lie on a single plane, infinitely many F matrices are consistent with
   the data.  The graffiti images used here are planar, so the F estimated by
   OpenCV is noisy — RANSAC compensates, but results are suboptimal.  In
   practice, these solvers are most useful when a reliable F is available from
   a scene with sufficient depth variation.

The workflow is:

1. Detect SIFT features and match with Lowe's ratio test.
2. Estimate the fundamental matrix F with OpenCV RANSAC.
3. Transform F to centred coordinates (the solver's native frame).
4. Run homlib LOMSAC with the 3PT solver.
5. Compare visually against the ground-truth homography.

"""

from time import time

import cv2
import numpy as np
import matplotlib.pyplot as plt

import homlib

###############################################################################
# Constants

IMAGE1_PATH = 'img/grafA.png'
IMAGE2_PATH = 'img/grafB.png'
GT_HOMOGRAPHY_PATH = 'img/graf_model.txt'

DESIRED_KEYPOINTS = 8000
SNN_THRESHOLD = 0.8
INLIER_THRESHOLD = 3.0

LOMSAC_KWARGS = dict(
    inlier_threshold=INLIER_THRESHOLD,
    final_least_squares=True,
    random_seed=2026,
    lo_starting_iterations=8,
    min_num_iterations=70,
    max_num_iterations=500,
)


###############################################################################
# Helpers

def center_points_and_get_transforms(points, width, height):
    """Shift points so that (0, 0) is at the image centre.

    Returns centred points, the centering matrix T, and its inverse T_inv.
    """
    shift = np.array([[width / 2.0], [height / 2.0]])
    points_centered = points - shift
    T = np.eye(3)
    T[0, 2] = -width / 2.0
    T[1, 2] = -height / 2.0
    T_inv = np.eye(3)
    T_inv[0, 2] = width / 2.0
    T_inv[1, 2] = height / 2.0
    return points_centered, T, T_inv


def denormalize_homography(H_centered, T1, T2_inv):
    """Convert a homography from centred back to pixel coordinates."""
    return T2_inv @ H_centered @ T1


def homography_error(H, H_gt):
    """Normalised Frobenius error between two homographies."""
    return np.linalg.norm(H / np.linalg.norm(H) - H_gt / np.linalg.norm(H_gt))


def decolorize(img):
    gray = cv2.cvtColor(img, cv2.COLOR_RGB2GRAY)
    return cv2.cvtColor(gray, cv2.COLOR_GRAY2RGB)


def draw_matches(kps1, kps2, tentatives, img1, img2,
                 H, H_gt, mask, title=""):
    if H is None:
        print("No homography found")
        return
    matches_mask = mask.ravel().tolist()
    h, w = img1.shape[:2]
    corners = np.float32(
        [[0, 0], [0, h - 1], [w - 1, h - 1], [w - 1, 0]]
    ).reshape(-1, 1, 2)
    dst_est = cv2.perspectiveTransform(corners, H)
    dst_gt = cv2.perspectiveTransform(corners, H_gt)
    img2_vis = decolorize(img2)
    img2_vis = cv2.polylines(
        img2_vis, [np.int32(dst_est)], True,
        (255, 0, 0), 3, cv2.LINE_AA)
    img2_vis = cv2.polylines(
        img2_vis, [np.int32(dst_gt)], True,
        (0, 255, 0), 3, cv2.LINE_AA)
    draw_params = dict(
        matchColor=(255, 255, 0), singlePointColor=None,
        matchesMask=matches_mask, flags=2)
    img_out = cv2.drawMatches(
        decolorize(img1), kps1, img2_vis, kps2,
        tentatives, None, **draw_params)
    plt.figure(figsize=(12, 8))
    plt.imshow(img_out)
    plt.title(title)
    plt.axis('off')


###############################################################################
# Load images and ground truth

img1 = cv2.cvtColor(cv2.imread(IMAGE1_PATH), cv2.COLOR_BGR2RGB)
img2 = cv2.cvtColor(cv2.imread(IMAGE2_PATH), cv2.COLOR_BGR2RGB)
H_gt = np.linalg.inv(np.loadtxt(GT_HOMOGRAPHY_PATH))

h1, w1 = img1.shape[:2]
h2, w2 = img2.shape[:2]

###############################################################################
# Step 1: Detect SIFT features and match

sift = cv2.SIFT_create(DESIRED_KEYPOINTS)
kps1, descs1 = sift.detectAndCompute(img1, None)
kps2, descs2 = sift.detectAndCompute(img2, None)

bf = cv2.BFMatcher()
knn_matches = bf.knnMatch(descs1, descs2, k=2)

tentatives = []
for m, n in knn_matches:
    if m.distance < SNN_THRESHOLD * n.distance:
        tentatives.append(m)

src_pts_px = np.float32(
    [kps1[m.queryIdx].pt for m in tentatives]).reshape(-1, 2)
dst_pts_px = np.float32(
    [kps2[m.trainIdx].pt for m in tentatives]).reshape(-1, 2)

print(f"{len(tentatives)} tentative matches after ratio test")

###############################################################################
# Step 2: Estimate the fundamental matrix F with OpenCV RANSAC

F_px, F_mask = cv2.findFundamentalMat(
    src_pts_px, dst_pts_px, cv2.FM_RANSAC,
    ransacReprojThreshold=1.0, confidence=0.999)
n_F_inliers = F_mask.sum()
print(f"Fundamental matrix: {n_F_inliers} / {len(tentatives)} inliers")

###############################################################################
# Step 3: Transform F to centred coordinates
#
# The solvers work in centred coordinates (origin at image centre).
# If T1, T2 are the centering matrices, then F_centered = T2^{-T} F_px T1^{-1}.

_, T1, T1_inv = center_points_and_get_transforms(src_pts_px.T, w1, h1)
_, T2, T2_inv = center_points_and_get_transforms(dst_pts_px.T, w2, h2)

F_centered = np.linalg.inv(T2).T @ F_px @ np.linalg.inv(T1)

###############################################################################
# Step 4: 3PT solver (3 points + F)
#
# This uses only point correspondences — no affine features needed.
# For the affine variant (1 AC + F), see :class:`~homlib.BarathVISAPP2016HAF`.

src_centered = (src_pts_px.T - np.array([[w1 / 2.0], [h1 / 2.0]]))
dst_centered = (dst_pts_px.T - np.array([[w2 / 2.0], [h2 / 2.0]]))

solver = homlib.BarathVISAPP20163PT(F=F_centered)
lomsac = homlib.LOMSAC(solver, **LOMSAC_KWARGS)

t = time()
estimate, stats = lomsac.estimate(src_centered, dst_centered)
elapsed = time() - t

H_3pt = denormalize_homography(estimate.homography, T1, T2_inv)

mask_3pt = np.array(
    [i in stats.inlier_indices for i in range(len(tentatives))],
    dtype=np.uint8)
print("\nhomlib 3PT + F:")
print(f"  H error     = {homography_error(H_3pt, H_gt):.4e}")
print(f"  Inliers     = {mask_3pt.sum()} / {len(tentatives)}")
print(f"  Time        = {elapsed:.3f} sec")

draw_matches(kps1, kps2, tentatives, img1, img2, H_3pt, H_gt, mask_3pt,
             title='homlib 3PT + F (blue = estimated, green = GT)')
