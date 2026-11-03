r"""
=============================
Real-data backend comparison
=============================

Compare RANSAC backends on real images with simulated radial distortion.

The pipeline:

1. Load two images and add a known radial distortion to the first.
2. Detect SIFT features, match with Lowe's ratio test.
3. Estimate a homography with OpenCV RANSAC (baseline) and every available
   HomLib backend, averaging over five runs per configuration.
4. Repeat with orientation-covariant features (SIFT angles).
5. (Optional) Repeat with affine features (AffNet + HardNet).
6. Print a summary table ranking all methods by accuracy and speed.

Backends: LOMSAC [1]_, SupeRANSAC [2]_, GC-RANSAC [3]_, MAGSAC++ [4]_.

Solvers: Wadenbäck et al. [5]_ (point-based), Valtonen Örnhag and
Adalbjörnsson [6]_ (orientation- and affine-covariant).

References
----------
.. [1] K. Lebeda, J. Matas, and O. Chum, "Fixing the Locally Optimized
   RANSAC," in Proc. BMVC, 2012.
.. [2] D. Barath, "SupeRANSAC: One RANSAC to Rule Them All," arXiv
   preprint arXiv:2506.04803, 2025.
.. [3] D. Barath and J. Matas, "Graph-Cut RANSAC," in Proc. CVPR, 2018.
.. [4] D. Barath, J. Noskova, M. Ivashechkin, and J. Matas, "MAGSAC++, a
   Fast, Reliable and Accurate Robust Estimator," in Proc. CVPR, 2020.
.. [5] M. Wadenbäck, M. Valtonen Örnhag, and J. Edstedt, "Radially Distorted
   Homographies, Revisited," in Proc. 3DV, 2026.
.. [6] M. Valtonen Örnhag and S. Adalbjörnsson, "Radial Distortion Homography
   Estimation from Affine-Covariant or Orientation-Covariant Features,"
   in Proc. ICPR, 2026.
"""

import random
import ctypes
from collections import defaultdict
from time import time

import cv2
import numpy as np
import matplotlib.pyplot as plt

import homlib

###############################################################################
# Configuration
# -------------

IMAGE1_PATH = "img/grafA.png"
IMAGE2_PATH = "img/grafB.png"
GT_HOMOGRAPHY_PATH = "img/graf_model.txt"

DESIRED_KEYPOINTS = 8000
DIST_COEFF_GT = -0.000003
SNN_THRESHOLD = 0.8
INLIER_THRESHOLD = 3.0
N_RUNS = 5

RANDOM_SEED = random.randint(0, ctypes.c_uint32(-1).value)

has_superansac = homlib._HAS_SUPERANSAC
has_gcransac = homlib._HAS_GCRANSAC
has_magsac = homlib._HAS_MAGSAC

print(
    f"Backends: LOMSAC=True, SupeRANSAC={has_superansac}, "
    f"GC-RANSAC={has_gcransac}, MAGSAC++={has_magsac}"
)

###############################################################################
# Image loading
# -------------
# We load two images and add a small radial distortion to the first one so
# that the distortion-aware solvers have something to estimate.


def load_images(img1_path, img2_path, dist_coeff):
    """Load images and apply radial distortion to the first."""
    img1 = cv2.cvtColor(cv2.imread(img1_path), cv2.COLOR_BGR2RGB)
    img2 = cv2.cvtColor(cv2.imread(img2_path), cv2.COLOR_BGR2RGB)

    K = np.eye(3)
    K[0, 2] = img2.shape[1] / 2.0
    K[1, 2] = img2.shape[0] / 2.0

    map_x, map_y = cv2.initUndistortRectifyMap(
        K, np.array([0, 0, 0, 0, 0, dist_coeff, 0, 0]),
        np.eye(3), K, (img2.shape[1], img2.shape[0]), cv2.CV_32FC1,
    )
    img1 = cv2.remap(img1, map_x, map_y, cv2.INTER_LINEAR)
    return img1, img2


img1, img2 = load_images(IMAGE1_PATH, IMAGE2_PATH, DIST_COEFF_GT)
H_gt = np.linalg.inv(np.loadtxt(GT_HOMOGRAPHY_PATH))
w1, h1 = img1.shape[1], img1.shape[0]
w2, h2 = img2.shape[1], img2.shape[0]

plt.figure()
plt.imshow(img1)
plt.title("Image 1 (distorted)")
plt.figure()
plt.imshow(img2)
plt.title("Image 2")


###############################################################################
# Feature detection
# -----------------


def detect_and_match_sift(img1, img2, max_kpts=DESIRED_KEYPOINTS):
    """Detect SIFT features and match with Lowe's ratio test."""
    sift = cv2.SIFT_create(max_kpts)
    kps1, descs1 = sift.detectAndCompute(img1, None)
    kps2, descs2 = sift.detectAndCompute(img2, None)

    bf = cv2.BFMatcher()
    knn_matches = bf.knnMatch(descs1, descs2, k=2)

    good = []
    ratios = []
    for m, n in knn_matches:
        if m.distance < SNN_THRESHOLD * n.distance:
            good.append(m)
            ratios.append(m.distance / n.distance)

    order = np.argsort(ratios)
    good = list(np.array(good)[order])
    return kps1, kps2, descs1, descs2, good


kps1, kps2, _, _, tentatives = detect_and_match_sift(img1, img2)
print(f"{len(tentatives)} tentative matches")


###############################################################################
# Helpers
# -------


def center_points(points, width, height):
    """Shift points so that (0, 0) is the image centre.  Returns (pts, T, T_inv)."""
    shift = np.array([[width / 2.0], [height / 2.0]])
    T = np.eye(3)
    T[0, 2] = -width / 2.0
    T[1, 2] = -height / 2.0
    T_inv = np.eye(3)
    T_inv[0, 2] = width / 2.0
    T_inv[1, 2] = height / 2.0
    return points - shift, T, T_inv


def homography_error(H, H_gt):
    """Frobenius error between normalised homographies."""
    return np.linalg.norm(H / np.linalg.norm(H) - H_gt / np.linalg.norm(H_gt))


def inlier_list(stats):
    """Extract inlier list from backend-specific statistics."""
    return getattr(stats, "inliers", None) or getattr(stats, "inlier_indices", [])


def make_backends(solver):
    """Create all available backends for *solver*."""
    opts = homlib.LORansacOptions()
    opts.squared_inlier_threshold = INLIER_THRESHOLD ** 2
    opts.final_least_squares = True
    opts.random_seed = RANDOM_SEED
    opts.lo_starting_iterations = 8
    opts.min_num_iterations = 70
    opts.max_num_iterations = 500

    backends = {"LOMSAC": homlib.LOMSAC(solver, options=opts)}
    if has_superansac:
        backends["SupeRANSAC"] = homlib.SupeRANSAC(solver, inlier_threshold=10.0)
    if has_gcransac:
        backends["GC-RANSAC"] = homlib.GCRANSAC(solver, threshold=INLIER_THRESHOLD)
    if has_magsac:
        backends["MAGSAC++"] = homlib.MAGSAC(solver, sigma_max=10.0)
    return backends


###############################################################################
# Point-based estimation
# ----------------------
# The one-sided solver assumes distortion on the *left* image, so we swap
# ``src`` and ``dst`` and invert the resulting homography.


def run_point(backend, src, dst):
    """Single point-based estimation run.  Returns (H, n_inliers, dist_err, elapsed)."""
    src_c, T1, _ = center_points(src, w1, h1)
    dst_c, _, T2_inv = center_points(dst, w2, h2)

    t0 = time()
    est, stats = backend.estimate(dst_c, src_c)
    elapsed = time() - t0

    H = T2_inv @ np.linalg.inv(est.homography) @ T1
    return H, len(inlier_list(stats)), est.distortion_parameter2, elapsed


src_pts = np.float32([kps1[m.queryIdx].pt for m in tentatives]).reshape(-1, 2).T
dst_pts = np.float32([kps2[m.trainIdx].pt for m in tentatives]).reshape(-1, 2).T

point_solver = homlib.Wadenback3DV2026OneSided()
point_backends = make_backends(point_solver)
point_results = defaultdict(lambda: defaultdict(list))

for name, backend in point_backends.items():
    for _ in range(N_RUNS):
        H, n_inl, d_err, elapsed = run_point(backend, src_pts, dst_pts)
        r = point_results[name]
        r["h_err"].append(homography_error(H, H_gt))
        r["d_err"].append(abs(d_err - DIST_COEFF_GT))
        r["inliers"].append(n_inl)
        r["time"].append(elapsed)
    r = point_results[name]
    print(
        f"  Point / {name:12s}  H err={np.median(r['h_err']):.2e}  "
        f"k err={np.median(r['d_err']):.2e}  "
        f"inliers={np.median(r['inliers']):.0f}  "
        f"time={np.median(r['time'])*1000:.1f} ms"
    )

###############################################################################
# Orientation-covariant estimation
# --------------------------------

src_ori = np.float32([kps1[m.queryIdx].angle for m in tentatives])
dst_ori = np.float32([kps2[m.trainIdx].angle for m in tentatives])
ori = np.vstack((src_ori, dst_ori)) / 180.0 * np.pi


def run_orientation(backend, src, dst, ori):
    """Single orientation-based estimation run."""
    src_c, T1, _ = center_points(src, w1, h1)
    dst_c, _, T2_inv = center_points(dst, w2, h2)

    t0 = time()
    est, stats = backend.estimate(src_c, dst_c, ori=ori)
    elapsed = time() - t0

    H = T2_inv @ est.homography @ T1
    return H, len(inlier_list(stats)), est.distortion_parameter, elapsed


ori_solver = homlib.ValtonenOrnhagICPR2026OrientationSingleSided()
ori_backends = make_backends(ori_solver)
ori_results = defaultdict(lambda: defaultdict(list))

for name, backend in ori_backends.items():
    for _ in range(N_RUNS):
        H, n_inl, d_err, elapsed = run_orientation(backend, src_pts, dst_pts, ori)
        r = ori_results[name]
        r["h_err"].append(homography_error(H, H_gt))
        r["d_err"].append(abs(d_err - DIST_COEFF_GT))
        r["inliers"].append(n_inl)
        r["time"].append(elapsed)
    r = ori_results[name]
    print(
        f"  Ori   / {name:12s}  H err={np.median(r['h_err']):.2e}  "
        f"k err={np.median(r['d_err']):.2e}  "
        f"inliers={np.median(r['inliers']):.0f}  "
        f"time={np.median(r['time'])*1000:.1f} ms"
    )


###############################################################################
# Affine-covariant estimation (optional)
# ---------------------------------------
# Requires ``kornia``, ``kornia_moons``, and ``torch``.

aff_results = defaultdict(lambda: defaultdict(list))

try:
    import kornia as K
    import kornia.feature as KF
    import torch
    from kornia_moons.feature import OpenCVDetectorWithAffNetKornia

    print("\nComputing AffNet + HardNet features ...")
    img1_t = K.image_to_tensor(img1, False).float() / 255.0
    img2_t = K.image_to_tensor(img2, False).float() / 255.0
    device = "cuda" if torch.cuda.is_available() else "cpu"
    img1_t, img2_t = img1_t.to(device), img2_t.to(device)

    detector = OpenCVDetectorWithAffNetKornia(
        cv2.SIFT_create(DESIRED_KEYPOINTS), max_kpts=DESIRED_KEYPOINTS
    )
    feature = KF.LocalFeature(detector, KF.LAFDescriptor(KF.HardNet(True)).eval())
    with torch.no_grad():
        lafs1, _, descs1 = feature(img1_t)
        lafs2, _, descs2 = feature(img2_t)

    lafs1 = np.squeeze(lafs1.cpu().numpy())
    descs1 = np.squeeze(descs1.cpu().numpy())
    lafs2 = np.squeeze(lafs2.cpu().numpy())
    descs2 = np.squeeze(descs2.cpu().numpy())

    knn = cv2.BFMatcher().knnMatch(descs1, descs2, k=2)
    tent_aff = [m for m, n in knn if m.distance < SNN_THRESHOLD * n.distance]

    xs = np.zeros((len(tent_aff), 2), dtype=np.float32)
    ys = np.zeros((len(tent_aff), 2), dtype=np.float32)
    A = np.zeros((len(tent_aff), 4), dtype=np.float32)
    for i, m in enumerate(tent_aff):
        L1, L2 = lafs1[m.queryIdx], lafs2[m.trainIdx]
        Amat = L2[:, :2] @ np.linalg.inv(L1[:, :2])
        xs[i] = [L1[0, 2], L1[1, 2]]
        ys[i] = [L2[0, 2], L2[1, 2]]
        A[i] = Amat.ravel()

    def run_affine(backend, xs, ys, A):
        src_c, T1, _ = center_points(xs.T, w1, h1)
        dst_c, _, T2_inv = center_points(ys.T, w2, h2)
        t0 = time()
        est, stats = backend.estimate(src_c, dst_c, A=A.T)
        elapsed = time() - t0
        H = T2_inv @ est.homography @ T1
        return H, len(inlier_list(stats)), est.distortion_parameter, elapsed

    aff_solver = homlib.ValtonenOrnhagICPR2026AffineSingleSided()
    aff_backends = make_backends(aff_solver)

    for name, backend in aff_backends.items():
        for _ in range(N_RUNS):
            H, n_inl, d_err, elapsed = run_affine(backend, xs, ys, A)
            r = aff_results[name]
            r["h_err"].append(homography_error(H, H_gt))
            r["d_err"].append(abs(d_err - DIST_COEFF_GT))
            r["inliers"].append(n_inl)
            r["time"].append(elapsed)
        r = aff_results[name]
        print(
            f"  Aff   / {name:12s}  H err={np.median(r['h_err']):.2e}  "
            f"k err={np.median(r['d_err']):.2e}  "
            f"inliers={np.median(r['inliers']):.0f}  "
            f"time={np.median(r['time'])*1000:.1f} ms"
        )

except ImportError:
    print("\nSkipping affine features (kornia/torch not available).")


###############################################################################
# OpenCV baseline
# ---------------

cv2_results = defaultdict(list)
src_cv = np.float32([kps1[m.queryIdx].pt for m in tentatives]).reshape(-1, 1, 2)
dst_cv = np.float32([kps2[m.trainIdx].pt for m in tentatives]).reshape(-1, 1, 2)

for _ in range(N_RUNS):
    t0 = time()
    H, mask = cv2.findHomography(src_cv, dst_cv, cv2.RANSAC, INLIER_THRESHOLD)
    elapsed = time() - t0
    cv2_results["h_err"].append(homography_error(H, H_gt))
    cv2_results["inliers"].append(int(mask.sum()))
    cv2_results["time"].append(elapsed)


###############################################################################
# Summary table
# -------------
# Median values over the five runs.  Sorted by homography error (best first).

rows = []

rows.append((
    "OpenCV RANSAC",
    np.median(cv2_results["h_err"]),
    float("nan"),
    np.median(cv2_results["inliers"]),
    np.median(cv2_results["time"]) * 1000,
))

for feat, results in [("Point", point_results), ("Ori", ori_results), ("Aff", aff_results)]:
    for name, r in results.items():
        rows.append((
            f"{feat} / {name}",
            np.median(r["h_err"]),
            np.median(r["d_err"]),
            np.median(r["inliers"]),
            np.median(r["time"]) * 1000,
        ))

rows.sort(key=lambda r: r[1])

print(f"\n{'Method':<30s} {'H err':>10s} {'k err':>10s} {'inliers':>8s} {'ms':>8s}")
print("-" * 70)
for method, h_err, d_err, inl, ms in rows:
    d_str = f"{d_err:10.2e}" if not np.isnan(d_err) else f"{'n/a':>10s}"
    print(f"{method:<30s} {h_err:10.2e} {d_str} {inl:8.0f} {ms:8.1f}")

best = rows[0]
fastest = min(rows, key=lambda r: r[4])
print(f"\nMost accurate: {best[0]} (H err = {best[1]:.2e})")
print(f"Fastest:       {fastest[0]} ({fastest[4]:.1f} ms)")


###############################################################################
# Visualization
# -------------
# Show inlier matches for the best-scoring configuration.


def decolorize(img):
    return cv2.cvtColor(cv2.cvtColor(img, cv2.COLOR_RGB2GRAY), cv2.COLOR_GRAY2RGB)


def draw_best(kps1, kps2, tentatives, img1, img2, H, H_gt, mask, title):
    if H is None:
        return
    h, w, _ = img1.shape
    corners = np.float32([[0, 0], [0, h - 1], [w - 1, h - 1], [w - 1, 0]]).reshape(-1, 1, 2)
    dst_est = cv2.perspectiveTransform(corners, H)
    dst_gt = cv2.perspectiveTransform(corners, H_gt)

    img2_vis = decolorize(img2)
    img2_vis = cv2.polylines(img2_vis, [np.int32(dst_est)], True, (255, 0, 0), 3, cv2.LINE_AA)
    img2_vis = cv2.polylines(img2_vis, [np.int32(dst_gt)], True, (0, 255, 0), 3, cv2.LINE_AA)

    draw_params = dict(matchColor=(255, 255, 0), singlePointColor=None,
                       matchesMask=mask.ravel().tolist(), flags=2)
    out = cv2.drawMatches(decolorize(img1), kps1, img2_vis, kps2, tentatives, None, **draw_params)
    plt.figure(figsize=(12, 8))
    plt.imshow(out)
    plt.title(title)
    plt.axis("off")


best_method = best[0]
print(f"\nVisualising best method: {best_method}")

if best_method.startswith("Point"):
    bname = best_method.split(" / ")[1]
    backend = point_backends[bname]
    H, n_inl, _, _ = run_point(backend, src_pts, dst_pts)
    mask = np.zeros(len(tentatives), dtype=np.uint8)
    draw_best(kps1, kps2, tentatives, img1, img2, H, H_gt, mask, best_method)
elif best_method.startswith("Ori"):
    bname = best_method.split(" / ")[1]
    backend = ori_backends[bname]
    H, n_inl, _, _ = run_orientation(backend, src_pts, dst_pts, ori)
    mask = np.zeros(len(tentatives), dtype=np.uint8)
    draw_best(kps1, kps2, tentatives, img1, img2, H, H_gt, mask, best_method)
elif best_method == "OpenCV RANSAC":
    H, mask = cv2.findHomography(src_cv, dst_cv, cv2.RANSAC, INLIER_THRESHOLD)
    draw_best(kps1, kps2, tentatives, img1, img2, H, H_gt, mask, best_method)

plt.show()
