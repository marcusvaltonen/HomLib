"""
Homography estimation with radial distortion
============================================

This example demonstrates:

1. Loading and pre‑processing two images (including adding a small radial distortion).
2. Detecting SIFT features and matching them with Lowe's ratio test.
3. Estimating a homography with OpenCV RANSAC and several homlib methods.
4. Comparing the results visually and numerically against ground truth.
5. Using AffNet + HardNet features and affine correspondences.

"""

import random
import ctypes
from time import time

import cv2
import numpy as np
import matplotlib.pyplot as plt

import homlib


###############################################################################
# Constants / Configuration

IMAGE1_PATH = 'img/grafA.png'
IMAGE2_PATH = 'img/grafB.png'
GT_HOMOGRAPHY_PATH = 'img/graf_model.txt'

DESIRED_KEYPOINTS = 8000          # Maximum number of SIFT keypoints
DIST_COEFF_GT = -0.000003         # Ground-truth radial distortion coefficient
SNN_THRESHOLD = 0.8               # Lowe's ratio‑test threshold
INLIER_THRESHOLD = 3.0            # Pixel reprojection error for RANSAC

# Create a reproducible but random seed
RANDOM_SEED = random.randint(0, ctypes.c_uint32(-1).value)

# homlib RANSAC options
RANSAC_OPTIONS = homlib.LORansacOptions()
RANSAC_OPTIONS.squared_inlier_threshold = INLIER_THRESHOLD ** 2
RANSAC_OPTIONS.final_least_squares = True
RANSAC_OPTIONS.random_seed = RANDOM_SEED
RANSAC_OPTIONS.lo_starting_iterations = 8
RANSAC_OPTIONS.min_num_iterations = 70
RANSAC_OPTIONS.max_num_iterations = 500


###############################################################################
# Image loading and pre‑processing

def load_images(img1_path, img2_path, dist_coeff):
    """Load the two images, apply radial distortion to the first one, and return RGB versions."""
    img1 = cv2.cvtColor(cv2.imread(img1_path), cv2.COLOR_BGR2RGB)
    img2 = cv2.cvtColor(cv2.imread(img2_path), cv2.COLOR_BGR2RGB)

    # Intrinsic matrix (assumed identity, principal point at image centre)
    K = np.eye(3)
    K[0, 2] = img2.shape[1] / 2.0   # cx = width/2
    K[1, 2] = img2.shape[0] / 2.0   # cy = height/2
    R = np.eye(3)

    # Apply radial distortion to img1 only (to simulate a real‑world effect)
    map_x, map_y = cv2.initUndistortRectifyMap(
        K,
        np.array([0, 0, 0, 0, 0, dist_coeff, 0, 0]),
        R,
        K,
        (img2.shape[1], img2.shape[0]),
        cv2.CV_32FC1,
    )
    img1 = cv2.remap(img1, map_x, map_y, cv2.INTER_LINEAR)

    return img1, img2


def load_ground_truth_homography(path):
    """Load ground‑truth homography from a text file and return its inverse."""
    return np.linalg.inv(np.loadtxt(path))


###############################################################################
# Let us look at the images we are working with

img1, img2 = load_images(IMAGE1_PATH, IMAGE2_PATH, DIST_COEFF_GT)
H_gt = load_ground_truth_homography(GT_HOMOGRAPHY_PATH)

# Display original images (optional)
plt.figure()
plt.imshow(img1)
plt.title('Image 1 (distorted)')
plt.figure()
plt.imshow(img2)
plt.title('Image 2')


###########################################################################
# Feature detection and matching


def detect_and_match_sift(img1, img2, max_kpts=DESIRED_KEYPOINTS):
    """
    Detect SIFT keypoints/descriptors and match them with a brute‑force
    matcher followed by Lowe's ratio test.
    """
    sift = cv2.SIFT_create(max_kpts)
    kps1, descs1 = sift.detectAndCompute(img1, None)
    kps2, descs2 = sift.detectAndCompute(img2, None)

    bf = cv2.BFMatcher()
    knn_matches = bf.knnMatch(descs1, descs2, k=2)

    # Apply ratio test and keep only good matches
    good_matches = []
    snn_ratios = []
    for m, n in knn_matches:
        if m.distance < SNN_THRESHOLD * n.distance:
            good_matches.append(m)
            snn_ratios.append(m.distance / n.distance)

    # Sort matches by their distance ratio (better matches first)
    sorted_indices = np.argsort(snn_ratios)
    good_matches = list(np.array(good_matches)[sorted_indices])

    return kps1, kps2, descs1, descs2, good_matches


#########################
# Visualization utilities

def decolorize(img):
    """Convert an RGB image to grayscale and back to RGB (for drawing)."""
    gray = cv2.cvtColor(img, cv2.COLOR_RGB2GRAY)
    return cv2.cvtColor(gray, cv2.COLOR_GRAY2RGB)


def draw_matches(kps1, kps2, tentatives, img1, img2, H, H_gt, mask, title=""):
    """
    Draw tentative correspondences and the estimated/ground‑truth homography
    quadrilaterals on the images.
    """
    if H is None:
        print("No homography found")
        return

    matches_mask = mask.ravel().tolist()

    h, w, _ = img1.shape
    corners = np.float32([[0, 0], [0, h - 1], [w - 1, h - 1], [w - 1, 0]]).reshape(-1, 1, 2)

    # Project image1 corners to image2 using estimated and ground‑truth homographies
    dst_est = cv2.perspectiveTransform(corners, H)
    dst_gt = cv2.perspectiveTransform(corners, H_gt)

    # Draw polygons on a copy of image2
    img2_vis = decolorize(img2)
    img2_vis = cv2.polylines(img2_vis, [np.int32(dst_est)], True, (255, 0, 0), 3, cv2.LINE_AA)  # blue
    img2_vis = cv2.polylines(img2_vis, [np.int32(dst_gt)], True, (0, 255, 0), 3, cv2.LINE_AA)   # green

    draw_params = dict(
        matchColor=(255, 255, 0),      # yellow
        singlePointColor=None,
        matchesMask=matches_mask,
        flags=2,
    )
    img_out = cv2.drawMatches(decolorize(img1), kps1, img2_vis, kps2, tentatives, None, **draw_params)

    plt.figure(figsize=(12, 8))
    plt.imshow(img_out)
    plt.title(title)
    plt.axis('off')


###############################
# Homography estimation helpers

def center_points_and_get_transforms(points, width, height):
    """
    Shift point coordinates so that (0,0) is at the image centre.
    Also return the matrices that transform original homogeneous coordinates
    to centered ones (T_center) and back to original (T_back).

    Parameters
    ----------
    points : np.ndarray, shape (2, N)
        Point coordinates (first row = x, second row = y).
    width, height : int
        Image dimensions.

    Returns
    -------
    points_centered : np.ndarray, shape (2, N)
        Shifted point coordinates.
    T_center : np.ndarray, shape (3, 3)
        Matrix mapping original homogeneous points to centered points.
    T_back : np.ndarray, shape (3, 3)
        Inverse of T_center, mapping centered points back to original.
    """
    shift = np.array([[width / 2.0], [height / 2.0]])
    points_centered = points - shift

    T_center = np.eye(3)
    T_center[0, 2] = -width / 2.0
    T_center[1, 2] = -height / 2.0

    T_back = np.eye(3)
    T_back[0, 2] = width / 2.0
    T_back[1, 2] = height / 2.0

    return points_centered, T_center, T_back


def denormalize_homography(H_centered, T1, T2_inv):
    """
    Convert a homography estimated in centered coordinates back to original image coordinates.
    """
    return T2_inv @ H_centered @ T1


def homography_error(H, H_gt):
    """Compute normalised error between two homographies."""
    H_norm = H / np.linalg.norm(H)
    H_gt_norm = H_gt / np.linalg.norm(H_gt)
    return np.linalg.norm(H_norm - H_gt_norm)


def print_estimation_results(method_name, distortion_parameter, H, H_gt):
    """Print the common information for every homography estimation method."""
    print(f"homlib {method_name}:")
    print(f"H error = {homography_error(H, H_gt):.4e}")
    print(f"Dist. coeff. error = {abs(distortion_parameter - DIST_COEFF_GT):.4e}")


#######################
# Let us first try OpenCV, which does not handle radial distortion.
# --- SIFT + BFMatcher + ratio test ---
kps1, kps2, _, _, tentatives = detect_and_match_sift(img1, img2)

def verify_cv2(kps1, kps2, tentatives, H_gt, inlier_thresh=INLIER_THRESHOLD):
    """Estimate homography with OpenCV RANSAC."""
    src_pts = np.float32([kps1[m.queryIdx].pt for m in tentatives]).reshape(-1, 1, 2)
    dst_pts = np.float32([kps2[m.trainIdx].pt for m in tentatives]).reshape(-1, 1, 2)

    H, mask = cv2.findHomography(src_pts, dst_pts, cv2.RANSAC, inlier_thresh)

    print(f"OpenCV RANSAC: {mask.astype(np.float32).sum():.0f} inliers")
    print(f"H error = {homography_error(H, H_gt):.4e}")
    return H, mask

t = time()
cv2_H, cv2_mask = verify_cv2(kps1, kps2, tentatives, H_gt)
print(f"{time() - t:.3f} sec (OpenCV RANSAC)\n")
draw_matches(kps1, kps2, tentatives, img1, img2, cv2_H, H_gt, cv2_mask,
             title='OpenCV RANSAC')

#######################
# The results are not terrible, but there is room for improvement. Let us take a look at
# the homlib versions.

def verify_homlib_point(kps1, kps2, tentatives, H_gt,
                        width1, height1, width2, height2,
                        options=RANSAC_OPTIONS):
    """Estimate homography with homlib using point correspondences (one‑sided radial distortion)."""
    src_pts = np.float32([kps1[m.queryIdx].pt for m in tentatives]).reshape(-1, 2).T
    dst_pts = np.float32([kps2[m.trainIdx].pt for m in tentatives]).reshape(-1, 2).T

    # Normalise coordinates to image centre for both images
    src_centered, T1, _ = center_points_and_get_transforms(src_pts, width1, height1)
    dst_centered, _, T2_inv = center_points_and_get_transforms(dst_pts, width2, height2)

    # This routine assumes the left images is distorted, hence switch places with src and dst
    estimate, stats = homlib.lomsac_nakano_icpr_2025_one_sided(
        dst_centered, src_centered, options
    )

    # Since we swapped dst and src, we seek the inverse
    H = denormalize_homography(np.linalg.inv(estimate.homography), T1, T2_inv)

    mask = np.array([i in stats.inlier_indices for i in range(len(tentatives))], dtype=np.uint8)

    print_estimation_results("point method", estimate.distortion_parameter2, H, H_gt)
    return H, mask


def verify_homlib_ori(kps1, kps2, tentatives, H_gt,
                      width1, height1, width2, height2,
                      options=RANSAC_OPTIONS):
    """Estimate homography with homlib using point correspondences + keypoint orientations."""
    src_pts = np.float32([kps1[m.queryIdx].pt for m in tentatives]).reshape(-1, 2).T
    dst_pts = np.float32([kps2[m.trainIdx].pt for m in tentatives]).reshape(-1, 2).T

    # Orientation data (in radians)
    src_ori = np.float32([kps1[m.queryIdx].angle for m in tentatives])
    dst_ori = np.float32([kps2[m.trainIdx].angle for m in tentatives])
    ori = np.vstack((src_ori, dst_ori)) / 180.0 * np.pi

    src_centered, T1, _ = center_points_and_get_transforms(src_pts, width1, height1)
    dst_centered, _, T2_inv = center_points_and_get_transforms(dst_pts, width2, height2)

    estimate, stats = homlib.lomsac_valtonenornhag_icpr_2026_one_sided_ori(
        src_centered, dst_centered, ori, options
    )

    H = denormalize_homography(estimate.homography, T1, T2_inv)

    mask = np.array([i in stats.inlier_indices for i in range(len(tentatives))], dtype=np.uint8)

    print_estimation_results("orientation method", estimate.distortion_parameter, H, H_gt)
    return H, mask


# --- 2. homlib: point‑based, one‑sided radial distortion ---
t = time()
homlib_H, homlib_mask = verify_homlib_point(
    kps1, kps2, tentatives, H_gt,
    img1.shape[1], img1.shape[0], img2.shape[1], img2.shape[0]
)
print(f"{time() - t:.3f} sec (homlib point)\n")
draw_matches(kps1, kps2, tentatives, img1, img2, homlib_H, H_gt, homlib_mask,
             title='homlib point')

# --- 3. homlib: point + orientation ---
t = time()
homlib_H_ori, homlib_mask_ori = verify_homlib_ori(
    kps1, kps2, tentatives, H_gt,
    img1.shape[1], img1.shape[0], img2.shape[1], img2.shape[0]
)
print(f"{time() - t:.3f} sec (homlib orientation)\n")
draw_matches(kps1, kps2, tentatives, img1, img2, homlib_H_ori, H_gt, homlib_mask_ori,
             title='homlib orientation')

############################################################
# We can use affine features too. Let's compute them first.


def verify_homlib_affine(src_pts, dst_pts, A, tentatives, H_gt,
                         width1, height1, width2, height2,
                         options=RANSAC_OPTIONS):
    """Estimate homography using affine correspondences (ACs)."""
    # src_pts and dst_pts are 2xN arrays
    src_centered, T1, _ = center_points_and_get_transforms(src_pts.T, width1, height1)
    dst_centered, _, T2_inv = center_points_and_get_transforms(dst_pts.T, width2, height2)

    estimate, stats = homlib.lomsac_valtonenornhag_icpr_2026_one_sided_affine(
        src_centered, dst_centered, A.T, options
    )

    H = denormalize_homography(estimate.homography, T1, T2_inv)

    mask = np.array([i in stats.inlier_indices for i in range(len(tentatives))], dtype=np.uint8)

    print_estimation_results("affine method", estimate.distortion_parameter, H, H_gt)
    return H, mask


# -----------------------------------------------------------------------------
# AffNet + HardNet feature handling
# -----------------------------------------------------------------------------

def compute_affnet_hardnet_features(img1, img2, desired_kpts):
    """
    Compute AffNet + HardNet features for the two input images.
    Returns local affine frames (LAFs) and descriptors for both images.
    """
    print("Computing AffNet + HardNet features (this may take a while) ...")
    import kornia as K
    import kornia.feature as KF
    import torch
    from kornia_moons.feature import OpenCVDetectorWithAffNetKornia

    # Convert images to torch tensors (already RGB, so no channel swap)
    img1_torch = K.image_to_tensor(img1, False).float() / 255.0
    img2_torch = K.image_to_tensor(img2, False).float() / 255.0

    device = "cuda" if torch.cuda.is_available() else "cpu"
    img1_torch = img1_torch.to(device)
    img2_torch = img2_torch.to(device)

    detector = OpenCVDetectorWithAffNetKornia(cv2.SIFT_create(desired_kpts), max_kpts=desired_kpts)
    descriptor = KF.LAFDescriptor(KF.HardNet(True)).eval()
    feature = KF.LocalFeature(detector, descriptor)

    with torch.no_grad():
        lafs1, _, descs1 = feature(img1_torch)
        lafs2, _, descs2 = feature(img2_torch)

    lafs1np = np.squeeze(lafs1.cpu().detach().numpy())
    descs1np = np.squeeze(descs1.cpu().detach().numpy())
    lafs2np = np.squeeze(lafs2.cpu().detach().numpy())
    descs2np = np.squeeze(descs2.cpu().detach().numpy())

    return lafs1np, descs1np, lafs2np, descs2np


def get_laf_centroids(lafs1, lafs2):
    """
    Extract point coordinates (centroids) from LAFs.
    Used only for visualization of matches.
    Returns two lists of (x, y) tuples.
    """
    kps1 = [(lafs1[i, 0, 2], lafs1[i, 1, 2]) for i in range(lafs1.shape[0])]
    kps2 = [(lafs2[i, 0, 2], lafs2[i, 1, 2]) for i in range(lafs2.shape[0])]
    return kps1, kps2


def get_affine_correspondences(lafs1, lafs2, tentatives):
    """
    Convert pairs of LAFs to affine correspondences (ACs).
    Returns:
        xs, ys : np.float32 arrays of shape (N, 2) containing centroids.
        A : np.float32 array of shape (N, 4) containing the affine transformation
            (flattened 2x2 matrix) between each pair of LAFs.
    """
    xs = np.zeros((len(tentatives), 2), dtype=np.float32)
    ys = np.zeros((len(tentatives), 2), dtype=np.float32)
    ACs = np.zeros((len(tentatives), 4), dtype=np.float32)
    for row, m in enumerate(tentatives):
        LAF1 = lafs1[m.queryIdx]
        LAF2 = lafs2[m.trainIdx]
        # Local affine transformation: A = LAF2 * inv(LAF1)
        A = np.matmul(LAF2[:, :2], np.linalg.inv(LAF1[:, :2]))
        xs[row, 0] = LAF1[0, 2]
        xs[row, 1] = LAF1[1, 2]
        ys[row, 0] = LAF2[0, 2]
        ys[row, 1] = LAF2[1, 2]
        ACs[row, 0] = A[0, 0]
        ACs[row, 1] = A[0, 1]
        ACs[row, 2] = A[1, 0]
        ACs[row, 3] = A[1, 1]
    return xs, ys, ACs

# --- 4. AffNet + HardNet features & affine correspondences ---
lafs1, descs1, lafs2, descs2 = compute_affnet_hardnet_features(
    img1, img2, DESIRED_KEYPOINTS
)

# Match descriptors with ratio test
bf = cv2.BFMatcher()
knn_matches = bf.knnMatch(descs1, descs2, k=2)
tentatives_aff = []
for m, n in knn_matches:
    if m.distance < SNN_THRESHOLD * n.distance:
        tentatives_aff.append(m)

# Convert LAFs to affine correspondences
xs, ys, A = get_affine_correspondences(lafs1, lafs2, tentatives_aff)

t = time()
homlib_aff_H, homlib_aff_mask = verify_homlib_affine(
    xs, ys, A, tentatives_aff, H_gt,
    img1.shape[1], img1.shape[0], img2.shape[1], img2.shape[0]
)
print(f"{time() - t:.3f} sec (homlib affine)\n")

# For visualization, create keypoints from LAF centroids
kps1_aff, kps2_aff = get_laf_centroids(lafs1, lafs2)
# Convert to cv2.KeyPoint objects for draw_matches
kps1_aff_cv = tuple(cv2.KeyPoint(x, y, 1) for x, y in kps1_aff)
kps2_aff_cv = tuple(cv2.KeyPoint(x, y, 1) for x, y in kps2_aff)

draw_matches(kps1_aff_cv, kps2_aff_cv, tentatives_aff,
             img1, img2, homlib_aff_H, H_gt, homlib_aff_mask,
             title='homlib affine')
