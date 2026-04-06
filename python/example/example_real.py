r"""
Real data experiment
====================

In this example, we show how to estimate a homography from images using
OpenCV's implementation of SIFT.
"""

import numpy as np
import matplotlib.pyplot as plt
import cv2
import homlib
from time import time


def decolorize(img):
    return  cv2.cvtColor(cv2.cvtColor(img,cv2.COLOR_RGB2GRAY), cv2.COLOR_GRAY2RGB)


def draw_matches(kps1, kps2, tentatives, img1, img2, H, mask):
    if H is None:
        print ("No homography found")
        return
    matchesMask = mask.ravel().tolist()
    h,w,ch = img1.shape
    pts = np.float32([ [0,0],[0,h-1],[w-1,h-1],[w-1,0] ]).reshape(-1,1,2)
    dst = cv2.perspectiveTransform(pts, H)
    #Ground truth transformation
    dst_GT = cv2.perspectiveTransform(pts, H_gt)
    img2_tr = cv2.polylines(decolorize(img2),[np.int32(dst)],True,(0,0,255),3, cv2.LINE_AA)
    img2_tr = cv2.polylines(deepcopy(img2_tr),[np.int32(dst_GT)],True,(0,255,0),3, cv2.LINE_AA)
    # Blue is estimated, green is ground truth homography
    draw_params = dict(matchColor = (255,255,0), # draw matches in yellow color
                   singlePointColor = None,
                   matchesMask = matchesMask, # draw only inliers
                   flags = 2)
    img_out = cv2.drawMatches(img1,kps1,img2_tr,kps2,tentatives,None,**draw_params)
    plt.figure(figsize = (12,8))
    plt.imshow(img_out)
    return


def verify_cv2(kps1, kps2, tentatives):
    src_pts = np.float32([ kps1[m.queryIdx].pt for m in tentatives ]).reshape(-1,1,2)
    dst_pts = np.float32([ kps2[m.trainIdx].pt for m in tentatives ]).reshape(-1,1,2)
    H, mask = cv2.findHomography(src_pts, dst_pts, cv2.RANSAC, 1.0)
    print (deepcopy(mask).astype(np.float32).sum(), 'inliers found')

  #  print(f"H error = {np.linalg.norm(H / np.linalg.norm(H) - H_gt / np.linalg.norm(H_gt)):4e}")

    return H, mask

def verify_homlib(kps1, kps2, tentatives, h1, w1, h2, w2):
    src_pts = np.float32([kps1[m.queryIdx].pt for m in tentatives]).reshape(-1, 2).T
    dst_pts = np.float32([kps2[m.trainIdx].pt for m in tentatives]).reshape(-1, 2).T

    # Shift s.t. (0,0) is in the center
    src_pts = src_pts - np.array([[h1, w1]]).T / 2.0
    dst_pts = dst_pts - np.array([[h2, w2]]).T / 2.0
    T1 = np.eye(3)
    T1[0, 2] = -h1 / 2.0
    T1[1, 2] = -w1 / 2.0
    T2inv = np.eye(3)
    T2inv[0, 2] = h2 / 2.0
    T2inv[1, 2] = w2 / 2.0

    options = homlib.LORansacOptions()
    options.squared_inlier_threshold = 5.0 ** 2
    options.min_num_iterations = 500
    options.final_least_squares = True
    estimate, stats = homlib.lomsac_kukelova_cvpr_2015_two_sided_equal(src_pts, dst_pts, options)

    H = T2inv @ estimate.homography @ T1

    print(estimate)
    mask = np.array([x in stats.inlier_indices for x in range(len(tentatives))], dtype=np.uint8)
    print(stats)
    print(estimate.distortion_parameter)
    print(estimate.distortion_parameter2)

    #print(f"H error = {np.linalg.norm(H / np.linalg.norm(H) - H_gt / np.linalg.norm(H_gt)):4e}")
    #print(f"Dist. coeff. error = {np.abs(estimate.distortion_parameter-dist_coeff_gt):4e}")

    return H, mask, estimate.distortion_parameter, stats.time

pose_data = np.load('../posedata.npy',allow_pickle=True).item()

from itertools import pairwise
from scipy.spatial.transform import Rotation



n_plane = np.array([0.5921775818439035, -0.3666382743010555, -0.7175667825220526])
n_norm = np.linalg.norm(n_plane)
d = 2.1878781959850406
n_plane /= n_norm
print(f'Plane normal: {n_plane}')
d /= n_norm

times = []
errors = []

for i, j in [(2512, 2513), (2513, 2514), (2514, 2515), (2517, 2518), (2518, 2519)]:  # pairwise(pose_data)
    print(f'Image pair: ({i}, {j})')
    p1 = pose_data[i]
    p2 = pose_data[j]

    R1 = Rotation.from_quat(p1['q']).as_matrix()
    R2 = Rotation.from_quat(p2['q']).as_matrix()
    t1 = p1['t']
    t2 = p2['t']
    K = p1['K']

    # Relative pose
    R = R2 @ R1.T
    t = t2 - R @ t1
    n = R1 @ n_plane

    H_gt = K @ (R + t @ n.T / d) @ np.linalg.inv(K)
    print('H_gt')
    print(H_gt/H_gt[2,2])
    img1 = cv2.cvtColor(cv2.imread(f'/home/elramav/Downloads/drive-download-20250816T102202Z-1-001/fisheye_grossmunster/images/DSC_{i}.JPG'), cv2.COLOR_BGR2RGB)
    img2 = cv2.cvtColor(cv2.imread(f'/home/elramav/Downloads/drive-download-20250816T102202Z-1-001/fisheye_grossmunster/images/DSC_{j}.JPG'), cv2.COLOR_BGR2RGB)

    #We will detect ORB features and match them with cross-check test
    det = cv2.SIFT_create(8000)
    kps1, descs1 = det.detectAndCompute(img1,None)
    kps2, descs2 = det.detectAndCompute(img2,None)

    bf = cv2.BFMatcher()

    SNN_threshold = 0.8
    matches = bf.knnMatch(descs1, descs2, k=2)

    # Apply ratio test
    snn_ratios = []
    tentatives = []
    for m, n in matches:
        if m.distance < SNN_threshold * n.distance:
            tentatives.append(m)
            snn_ratios.append(m.distance / n.distance)

    sorted_indices = np.argsort(snn_ratios)
    tentatives = list(np.array(tentatives)[sorted_indices])

    print(f'Tentative: {len(tentatives)}')

    #Now, some visualization from OpenCV tutorial
    #https://docs.opencv.org/3.0-beta/doc/py_tutorials/py_feature2d/py_feature_homography/py_feature_homography.html
    #We will draw correspondences found and the geometric transformation between the images.
    from copy import deepcopy

    H, homlib_mask, kappa, exec_time = verify_homlib(kps1, kps2, tentatives, img1.shape[1], img1.shape[0], img2.shape[1], img2.shape[0])
    print(exec_time / 1e6, ' ms homlib')
    print(H / H[2,2])

    print(f"H error = {np.linalg.norm(H / np.linalg.norm(H) - H_gt / np.linalg.norm(H_gt)):4e}")


    times.append(exec_time / 1e6)
    errors.append(np.linalg.norm(H / np.linalg.norm(H) - H_gt / np.linalg.norm(H_gt)))

    #draw_matches(kps1, kps2, tentatives, img1, img2, H, homlib_mask)

    #plt.show()

print(f'time: {np.mean(times)}, (median {np.median(times)})')
print(f'error: {np.mean(errors)}, (median {np.median(errors)})')
