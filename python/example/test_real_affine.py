import numpy as np
import matplotlib.pyplot as plt
import cv2
import pygcransac
from copy import deepcopy
import kornia as K
import kornia.feature as KF
import torch
from kornia_moons.feature import OpenCVDetectorWithAffNetKornia
from kornia_moons.viz import draw_LAF_matches, visualize_LAF
from time import time

import homlib

from enum import Enum


class Solver(Enum):
    PointBased = 0
    SIFTBased = 1
    AffineBased = 2


class Sampler(Enum):
    Uniform = 0
    PROSAC = 1
    PNAPSAC = 2
    NGRANSAC = 3
    ARSampler = 4


# Function to load the image into a pytorch tensor
def load_torch_image(fname):
    img = K.image_to_tensor(cv2.imread(fname), False).float() / 255.
    img = K.color.bgr_to_rgb(img)
    return img


# Deciding about the device used. Prefer CUDA if available.
device = "cpu"
if torch.cuda.is_available():
    device = "cuda"

# The number of keypoints to be detected
desired_kpts = 2000

# Loading the images
print("Loading images")

# img1 = load_torch_image('img/grafA.png')
img1 = load_torch_image('/home/elramav/Downloads/drive-download-20250816T102202Z-1-001/fisheye_grossmunster/images/DSC_2515.JPG')
img2 = load_torch_image('/home/elramav/Downloads/drive-download-20250816T102202Z-1-001/fisheye_grossmunster/images/DSC_2517.JPG')


H_gt = np.linalg.inv(np.loadtxt('img/graf_model.txt'))
print(H_gt / H_gt[2, 2])

# Sending the images to the device
img1.to(device)
img2.to(device)

# Detecting DoG + HardNet + AffNet features
print("Computing features")
kornia_cv2dogaffnet = OpenCVDetectorWithAffNetKornia(cv2.SIFT_create(desired_kpts), max_kpts=desired_kpts)
dogaffnethardnet = KF.LocalFeature(kornia_cv2dogaffnet, KF.LAFDescriptor(KF.HardNet(True))).eval()

# Detecting features in the source image
lafs1, r1, descs1 = dogaffnethardnet(img1)
# Detecting features in the destination image
lafs2, r2, descs2 = dogaffnethardnet(img2)

# Visualizing the found affine correspondences
print("Visualizing")
visualize_LAF(img1, lafs1, 0, 'y', figsize=(8,6))
visualize_LAF(img2, lafs2, 0, 'y', figsize=(8,6))
plt.show()
# Initialize a brute-force matcher to match the descriptors
print("matching")
bf = cv2.BFMatcher()

# Send the detected keypoints and other variables to the CPU and convert them to numpy array
descs1np = np.squeeze(descs1.cpu().detach().numpy())
descs2np = np.squeeze(descs2.cpu().detach().numpy())
lafs1np = np.squeeze(lafs1.cpu().detach().numpy())
lafs2np = np.squeeze(lafs2.cpu().detach().numpy())
r1np = np.squeeze(r1.cpu().detach().numpy())
r2np = np.squeeze(r2.cpu().detach().numpy())

# The threshold for the SNN ratio test
SNN_threshold = 0.8

# Applying brute-force matcher
matches = bf.knnMatch(descs1np, descs2np, k=2)

# Apply the SNN ratio test
snn_ratios = []
tentatives = []
for m, n in matches:
    if m.distance < SNN_threshold * n.distance:
        tentatives.append(m)
        snn_ratios.append(m.distance / n.distance)

# Sort the keypoints based on the SNN ratio.
# This is used in many samplers, e.g., PROSAC, NG-RANSAC's, AR-Sampler
sorted_indices = np.argsort(snn_ratios)
tentatives = list(np.array(tentatives)[sorted_indices])

print(f"{len(tentatives)} tentative correspondences are found.")


# A function to convert to local affine frames (LAFs) to their centroids to obtain simple keypoints
def get_coordinates(lafs1, lafs2):
    kps1 = [[lafs1[i, 0, 2], lafs1[i, 1, 2]] for i in range(lafs1.shape[0])]
    kps2 = [[lafs2[i, 0, 2], lafs2[i, 1, 2]] for i in range(lafs2.shape[0])]
    return kps1, kps2


# A function to run OpenCV's RANSAC on the point correspondences extracted from the affine frames
def verify_cv2_homography(kps1, kps2, tentatives):
    # Copy the coordinates in the source image selected by the tentative correspondences
    src_pts = np.float32([kps1[m.queryIdx] for m in tentatives]).reshape(-1, 1, 2)
    # Copy the coordinates in the destination image selected by the tentative correspondences
    dst_pts = np.float32([kps2[m.trainIdx] for m in tentatives]).reshape(-1, 1, 2)
    # Apply OpenCV's RANSAC
    H, mask = cv2.findHomography(src_pts, dst_pts, cv2.RANSAC, 1.5)
    print(deepcopy(mask).astype(np.float32).sum(), 'inliers found')
    return H, mask


def verify_homlib_homography(src_pts, dst_pts, h1, w1, h2, w2):
    # Run LO-RANSAC with the solver
    options = homlib.LORansacOptions()
    options.squared_inlier_threshold = 10.0 ** 2
    options.final_least_squares = True

    # Shift s.t. (0,0) is in the center
    src_pts = np.array(src_pts).T - np.array([[h1, w1]]).T / 2.0
    dst_pts = np.array(dst_pts).T - np.array([[h2, w2]]).T / 2.0
    T1 = np.eye(3)
    T1[0, 2] = -h1 / 2.0
    T1[1, 2] = -w1 / 2.0
    T2inv = np.eye(3)
    T2inv[0, 2] = h2 / 2.0
    T2inv[1, 2] = w2 / 2.0

    estimate, stats = homlib.lomsac_wadenback_2025_two_sided_equal(src_pts, dst_pts, options)
    print(stats)
    print(estimate)

    mask = np.zeros(x.shape[0], dtype=bool)
    mask[stats.inlier_indices] = 1
    print(deepcopy(mask).astype(np.float32).sum(), 'inliers found')

    H = T2inv @ estimate.homography @ T1

    return H, mask


# Extracting the point correspondences to OpenCV RANSAC can run
kps1, kps2 = get_coordinates(lafs1np, lafs2np)

print(f"H gt =\n{H_gt / H_gt[2, 2]}\n")

# The time before starting the estimation with OpenCV's RANSAC
t = time()
# Running OpenCV's RANSAC for fundamental matrix estimation
cv2_H, cv2_mask = verify_cv2_homography(kps1, kps2, tentatives)
# Measuring the run-time
print(time() - t, 'sec cv2')


# The time before starting the estimation with OpenCV's RANSAC
t = time()
gc_H, gc_mask = verify_homlib_homography(kps1, kps2,
                                                    # The matches containing the indices of the matched keypoints
                                                    img1.shape[3],  # The width of the source image
                                                    img1.shape[2],  # The height of the source image
                                                    img2.shape[3],  # The width of the destination image
                                                    img2.shape[2],  # The height of the destination image
)
# The id of the used sampler. 0 - uniform, 1 - PROSAC, 3 - NG-RANSAC's sampler, 4 - AR-Sampler
# Measuring the run-time
print(time() - t, 'sec pygcransac-ac')
print(f"H pygcransac-ac =\n{gc_H / gc_H[2, 2]}\n")

tent_idxs = torch.from_numpy(np.array([[m.queryIdx, m.trainIdx] for m in tentatives]))

inlier_color = (0.2, 1, 0.2)
inlier_color = None
with torch.no_grad():
    draw_LAF_matches(lafs1, lafs2, tent_idxs,
                     img1, img2,
                     inlier_mask=cv2_mask,
                     draw_dict={"inlier_color": inlier_color,
                                "tentative_color": None,
                                "feature_color": None,
                                "vertical": False}, H=cv2_H)

    draw_LAF_matches(lafs1, lafs2, tent_idxs,
                     img1, img2,
                     inlier_mask=gc_mask,
                     draw_dict={"inlier_color": inlier_color,
                                "tentative_color": None,
                                "feature_color": None,
                                "vertical": False}, H=gc_H)

    draw_LAF_matches(lafs1, lafs2, tent_idxs,
                     img1, img2,
                     inlier_mask=homlib_mask,
                     draw_dict={"inlier_color": inlier_color,
                                "tentative_color": None,
                                "feature_color": None,
                                "vertical": False}, H=homlib_H)

plt.show()

