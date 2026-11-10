"""
Ground truth and matches of the Grossmünster Church pairs, for grossmunster_benchmark.py.

Wadenbäck et al. (3DV 2026, Section 5) evaluate the two-sided solvers on six pairs of the fisheye
images of the Grossmünster Church (RadialSfM, Larsson et al., ECCV 2020), with a ground-truth
homography from a 3D reconstruction. The script rebuilds that ground truth:

1. A 3D reconstruction (pycolmap, incremental SfM) of the images DSC_2510 to DSC_2522, with one
   OPENCV_FISHEYE camera shared by all the images (fx, fy, cx, cy and 4 coefficients of the distorted
   angle), from COLMAP's SIFT features at full resolution, exhaustively matched. By default
   (``--intrinsics paper``), the camera is the one the paper's authors gave (issue #42 of HomLib),
   fixed: f = 1917.57, (cx, cy) = (2144, 1424), the image center, and k1 = -0.04162, k2 = -0.00726
   (their "SIMPLE_RADIAL" coefficients fit OPENCV_FISHEYE: with them, COLMAP's RADIAL model
   registers only 12 of the 13 images and 1236 points, OPENCV_FISHEYE all of them and 10980 points).
   ``--intrinsics refined`` refines it in the reconstruction instead, which halves the reprojection
   errors (mean 0.48 px, against 0.88 px), mostly by moving the principal point by (13, 20) px.
   The script prints the reprojection errors. The reconstruction is saved in ``--sfm_dir`` and reused
   by later runs (whatever their ``--images``, ``--intrinsics`` and ``--seed``: give them another
   ``--sfm_dir``).
2. The plane of the wall with the two towers, fitted by RANSAC on the 3D points, then by least
   squares on its inliers. The wall isn't flat (pilasters, windows, cornices): the inlier threshold,
   a small share of the median depth of the points (``--plane_threshold``), keeps its flat parts.
3. The ground-truth homography of each pair, as in the paper: with the poses x_i = R_i X + t_i, the
   relative pose R = R2 R1^T, t = t2 - R t1, and the plane n'^T x = d in the coordinates of camera 1,
   H = K (R + t n'^T / d) K^-1 (the paper writes the plane n'^T x + d = 0, hence its minus sign). K
   is the camera's focal lengths: H maps the undistorted points, centered on the principal point,
   in pixels of a pinhole camera with the focal lengths of the fisheye camera at its center. These
   are the coordinates of HomLib's division model, whose undistortion keeps the scale at the
   distortion center: the solvers estimate H in them, from the points centered on the principal
   point.
4. SIFT matches of each pair (OpenCV SIFT, 8000 keypoints, Lowe's ratio test at 0.8, as the paper;
   hpatches_features.py), saved with the ground truth in a ``.npz`` file.

The script checks the ground truth on two sets of correspondences, mapped with it and with the
calibration (the fisheye undistortion of image 1, H, the fisheye distortion of image 2), in pixels of
image 2:
- the tracks of the reconstruction on the plane seen in both images: their residuals, against those
  of the homography fitted to them by least squares, which is the best that a homography can do;
  close values show that the poses and the calibration aren't what limits the ground truth, but the
  flatness of the wall;
- the SIFT matches: their number, and those within 1 and 3 px.
It also prints these numbers for the other pairs of consecutive images, which shows how the pairs
were chosen (more than 100 tentative matches; see PAIRS in grossmunster_data.py), and saves those
with more than 100 matches too, and the other candidate for the paper's sixth pair (OTHER_PAIRS), which
``grossmunster_benchmark.py --pairs`` can select.

Pixel coordinates: OpenCV puts (0, 0) at the center of the top-left pixel, COLMAP at the top-left
corner of the image. The saved points are in OpenCV's convention, centered on the principal point.

Data
----
The images of RadialSfM's "Grossmunster Church" (https://github.com/vlarsson/radialsfm, Google
Drive link in its README; Google Drive refuses scripted downloads), unzipped: ``--data_dir`` is the
directory that contains ``images/``.

Requirements
------------
pip install pycolmap opencv-python tqdm

Usage
-----
    python grossmunster_ground_truth.py --data_dir data/fisheye_grossmunster --sfm_dir data/grossmunster_sfm \
        --output data/grossmunster_sift.npz
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pycolmap

from grossmunster_data import OTHER_PAIRS, PAIRS, image_name, pair_key, save_matches

# hpatches_features.py's SIFT matcher
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "hpatches"))
from hpatches_features import sift_matcher  # noqa: E402

# The images of the reconstruction, DSC_2510 to DSC_2522
IMAGES = (2510, 2522)
# The paper keeps the pairs with more than 100 tentative matches
MIN_MATCHES = 100


# The camera of the paper (issue #42 of HomLib), as OPENCV_FISHEYE parameters (fx, fy, cx, cy, k1..k4),
# in COLMAP's convention (the origin at the top-left corner of the image)
PAPER_CAMERA = [1917.57, 1917.57, 2144.0, 1424.0, -0.04162, -0.00726, 0.0, 0.0]


def reconstruct(image_dir, sfm_dir, names, seed, intrinsics):
    """The reconstruction of the images, with one OPENCV_FISHEYE camera, the paper's (fixed) or refined
    (intrinsics); reused if sfm_dir has one."""
    model_dir = sfm_dir / "0"
    if (model_dir / "images.bin").exists():
        print(f"Reusing the reconstruction in {model_dir}")
        rec = pycolmap.Reconstruction(model_dir)
        fixed = np.allclose(next(iter(rec.cameras.values())).params, PAPER_CAMERA)
        if fixed != (intrinsics == "paper"):
            sys.exit(f"{model_dir} has a camera {'fixed to' if fixed else 'other than'} the paper's: give "
                     f"--intrinsics {'paper' if fixed else 'refined'}, or another --sfm_dir")
        return rec
    sfm_dir.mkdir(parents=True, exist_ok=True)
    database = sfm_dir / "database.db"
    database.unlink(missing_ok=True)
    pycolmap.set_random_seed(seed)
    extraction = pycolmap.FeatureExtractionOptions()
    extraction.max_image_size = 100000   # full resolution
    extraction.sift.max_num_features = 16384
    reader = pycolmap.ImageReaderOptions(camera_model="OPENCV_FISHEYE")
    if intrinsics == "paper":
        reader.camera_params = ",".join(map(str, PAPER_CAMERA))
    pycolmap.extract_features(database, image_dir, image_names=names, camera_mode=pycolmap.CameraMode.SINGLE,
                              reader_options=reader, extraction_options=extraction)
    pycolmap.match_exhaustive(database)
    options = pycolmap.IncrementalPipelineOptions()
    if intrinsics == "paper":
        options.ba_refine_focal_length = options.ba_refine_principal_point = options.ba_refine_extra_params = False
    else:
        options.ba_refine_principal_point = True
    reconstructions = pycolmap.incremental_mapping(database, image_dir, sfm_dir, options=options)
    best = max(reconstructions.values(), key=lambda r: r.num_reg_images())
    if best.num_reg_images() < len(names):
        print(f"Warning: the reconstruction has {best.num_reg_images()} of the {len(names)} images")
    best.write(model_dir)
    return best


def reprojection_errors(rec):
    """The reprojection errors of every observation, in pixels, per image name"""
    cam = next(iter(rec.cameras.values()))
    errors = {}
    for im in rec.images.values():
        obs = [(p.xy, rec.points3D[p.point3D_id].xyz) for p in im.points2D if p.has_point3D()]
        xy, X = (np.array(v) for v in zip(*obs))
        T = im.cam_from_world()
        Xc = X @ T.rotation.matrix().T + T.translation
        errors[im.name] = np.linalg.norm(np.array(cam.img_from_cam(Xc)) - xy, axis=1)
    return errors


def fit_plane(X, threshold, iterations=20000, seed=0):
    """(n, d, inliers): the plane n^T X = d (|n| = 1) with the most points within the threshold
    (RANSAC), refined by least squares on its inliers."""
    rng = np.random.default_rng(seed)
    best, n, d = -1, None, None
    for _ in range(iterations):
        s = X[rng.choice(len(X), 3, replace=False)]
        m = np.cross(s[1] - s[0], s[2] - s[0])
        if np.linalg.norm(m) < 1e-12:
            continue
        m /= np.linalg.norm(m)
        count = np.sum(np.abs(X @ m - m @ s[0]) < threshold)
        if count > best:
            best, n, d = count, m, m @ s[0]
    for _ in range(10):
        inliers = np.abs(X @ n - d) < threshold
        c = X[inliers].mean(0)
        n = np.linalg.svd(X[inliers] - c)[2][2]
        d = n @ c
    return n, d, np.abs(X @ n - d) < threshold


class Camera:
    """The calibration, in OpenCV's pixel convention: centered points <-> pixels, and the pinhole
    coordinates (undistorted, centered, in pixels of focal lengths fx, fy) of the fisheye pixels."""

    def __init__(self, cam):
        self.cam = cam
        self.fx, self.fy, cx, cy = cam.params[:4]
        self.center = np.array([cx - 0.5, cy - 0.5])
        self.K = np.diag([self.fx, self.fy, 1.0])

    def undistort(self, p):
        """Pinhole coordinates of the pixels p (Nx2)"""
        return np.array(self.cam.cam_from_img(p + 0.5)).reshape(-1, 2) * [self.fx, self.fy]

    def distort(self, u):
        """Pixels of the pinhole coordinates u (Nx2)"""
        rays = np.c_[u / [self.fx, self.fy], np.ones(len(u))]
        return np.array(self.cam.img_from_cam(rays)).reshape(-1, 2) - 0.5

    def map(self, H, p):
        """The pixels p of image 1 mapped to image 2 with the homography H of the pinhole
        coordinates; NaN for those mapped behind the camera"""
        q = H @ np.vstack([self.undistort(p).T, np.ones(len(p))])
        ok = q[2] > 0
        out = np.full((len(p), 2), np.nan)
        out[ok] = self.distort((q[:2, ok] / q[2, ok]).T)
        return out


def ground_truth(rec, name1, name2, n, d, K):
    """H = K (R + t n'^T / d') K^-1, from the poses of the images and the plane n^T X = d"""
    T1 = rec.find_image_with_name(name1).cam_from_world()
    T2 = rec.find_image_with_name(name2).cam_from_world()
    R1, t1, R2, t2 = T1.rotation.matrix(), T1.translation, T2.rotation.matrix(), T2.translation
    R = R2 @ R1.T
    t = t2 - R @ t1
    n1 = R1 @ n               # the plane in camera 1: n1^T x = d + n1^T t1
    return K @ (R + np.outer(t, n1) / (d + n1 @ t1)) @ np.linalg.inv(K)


def plane_tracks(rec, name1, name2, inliers):
    """The pixels (OpenCV's convention) of the tracks of the plane points seen in both images"""
    im1, im2 = rec.find_image_with_name(name1), rec.find_image_with_name(name2)
    ids = rec.point3D_ids()
    on_plane = {i for i, ok in zip(ids, inliers) if ok}
    xy2 = {p.point3D_id: p.xy for p in im2.points2D if p.has_point3D() and p.point3D_id in on_plane}
    pairs = [(p.xy, xy2[p.point3D_id]) for p in im1.points2D if p.has_point3D() and p.point3D_id in xy2]
    p1, p2 = (np.array(v).reshape(-1, 2) - 0.5 for v in zip(*pairs)) if pairs else (np.zeros((0, 2)),) * 2
    return p1, p2


def fitted_homography(u1, u2):
    """The homography of the pinhole coordinates fitted to u1 -> u2 by linear least squares (DLT,
    normalized), the best that a homography can do on these correspondences"""
    def normalize(u):
        c, s = u.mean(0), np.sqrt(2) / np.mean(np.linalg.norm(u - u.mean(0), axis=1))
        return np.array([[s, 0, -s * c[0]], [0, s, -s * c[1]], [0, 0, 1]])
    T1, T2 = normalize(u1), normalize(u2)
    a = (T1 @ np.c_[u1, np.ones(len(u1))].T).T
    b = (T2 @ np.c_[u2, np.ones(len(u2))].T).T
    rows = [np.r_[np.zeros(3), -b[i, 2] * a[i], b[i, 1] * a[i]] for i in range(len(a))]
    rows += [np.r_[b[i, 2] * a[i], np.zeros(3), -b[i, 0] * a[i]] for i in range(len(a))]
    H = np.linalg.inv(T2) @ np.linalg.svd(np.array(rows))[2][-1].reshape(3, 3) @ T1
    # The sign that maps the points in front of the camera
    return H * np.sign(np.mean(H[2] @ np.c_[u1, np.ones(len(u1))].T))


def main():
    parser = argparse.ArgumentParser(description="Ground truth and matches of the Grossmünster Church pairs")
    parser.add_argument("--data_dir", required=True, help="RadialSfM's fisheye_grossmunster (with images/)")
    parser.add_argument("--sfm_dir", required=True,
                        help="Directory of the reconstruction: created, or reused (whatever --images and --seed)")
    parser.add_argument("--output", required=True, help="The matches file to write (.npz)")
    parser.add_argument("--plane_threshold", type=float, default=0.003,
                        help="Inlier threshold of the plane, as a share of the median depth of the 3D points "
                             "in all the cameras (default: 0.003)")
    parser.add_argument("--num_features", type=int, default=8000,
                        help="SIFT keypoints per image (default: 8000, as the paper)")
    parser.add_argument("--ratio", type=float, default=0.8, help="Ratio test threshold (default: 0.8)")
    parser.add_argument("--images", type=int, nargs=2, default=IMAGES, metavar=("FIRST", "LAST"),
                        help="The images of the reconstruction, DSC_<FIRST> to DSC_<LAST> (default: %(default)s)")
    parser.add_argument("--seed", type=int, default=1, help="Random seed of the reconstruction and of the plane")
    parser.add_argument("--intrinsics", choices=["paper", "refined"], default="paper",
                        help="The camera: the paper's, fixed (default), or refined in the reconstruction")
    args = parser.parse_args()
    image_dir = Path(args.data_dir) / "images"
    numbers = range(args.images[0], args.images[1] + 1)

    # 1. The reconstruction
    rec = reconstruct(image_dir, Path(args.sfm_dir), [image_name(i) for i in numbers], args.seed, args.intrinsics)
    camera = Camera(next(iter(rec.cameras.values())))
    errors = reprojection_errors(rec)
    all_errors = np.concatenate(list(errors.values()))
    print(f"\nReconstruction: {rec.num_reg_images()} images, {rec.num_points3D()} points, "
          f"{len(all_errors)} observations; camera {camera.cam.model.name} "
          f"{np.array2string(camera.cam.params, precision=6, max_line_width=200)}")
    print(f"Reprojection error (px): mean {all_errors.mean():.3f}, median {np.median(all_errors):.3f}, "
          f"90th percentile {np.percentile(all_errors, 90):.3f}")
    for name, e in sorted(errors.items()):
        print(f"  {name}: {len(e):5d} observations, mean {e.mean():.3f}, median {np.median(e):.3f}")

    # 2. The plane
    X = np.array([rec.points3D[i].xyz for i in rec.point3D_ids()])
    depths = np.concatenate([X @ im.cam_from_world().rotation.matrix()[2] + im.cam_from_world().translation[2]
                             for im in rec.images.values()])
    threshold = args.plane_threshold * np.median(depths)
    n, d, inliers = fit_plane(X, threshold, seed=args.seed)
    print(f"\nPlane: {inliers.sum()} of {len(X)} points within {threshold:.4f} (RMS distance "
          f"{np.sqrt(np.mean((X[inliers] @ n - d) ** 2)):.4f}; median depth {np.median(depths):.3f})")

    # 3 and 4. The ground truth and the matches, and their checks
    match = sift_matcher(args.num_features, args.ratio)
    candidates = sorted(set(PAIRS) | set(OTHER_PAIRS) | {(i, i + 1) for i in numbers[:-1]})
    print(f"\n{'Pair':<10s} {'Plane tracks':>12s} {'GT median/p90':>14s} {'LS fit median/p90':>18s}   "
          f"{'SIFT matches':>12s} {'< 1 px':>7s} {'< 3 px':>7s}")
    pairs = {}
    for a, b in candidates:
        name1, name2 = image_name(a), image_name(b)
        H = ground_truth(rec, name1, name2, n, d, camera.K)
        t1, t2 = plane_tracks(rec, name1, name2, inliers)
        if len(t1) >= 8:
            r_gt = np.linalg.norm(camera.map(H, t1) - t2, axis=1)
            H_fit = fitted_homography(camera.undistort(t1), camera.undistort(t2))
            r_fit = np.linalg.norm(camera.map(H_fit, t1) - t2, axis=1)
        else:   # too few tracks to check the ground truth (or a reconstruction in pieces)
            r_gt = r_fit = np.full(1, np.nan)
        p1, p2, scores = match(str(image_dir / name1), str(image_dir / name2))
        r = np.linalg.norm(camera.map(H, p1) - p2, axis=1)
        mark = " *" if (a, b) in PAIRS else ""
        print(f"{pair_key(a, b):<10s} {len(t1):12d} {np.median(r_gt):6.2f} {np.percentile(r_gt, 90):6.2f} "
              f"{np.median(r_fit):9.2f} {np.percentile(r_fit, 90):6.2f}     {len(p1):12d} {np.sum(r < 1):7d} "
              f"{np.sum(r < 3):7d}{mark}")
        if (a, b) in PAIRS or (a, b) in OTHER_PAIRS or len(p1) > MIN_MATCHES:
            pairs[pair_key(a, b)] = dict(x=(p1 - camera.center).T, y=(p2 - camera.center).T, H=H, scores=scores)
    print(f"(* the benchmark's default pairs; residuals in pixels of the second image. The file has these pairs, "
          f"{', '.join(pair_key(*p) for p in OTHER_PAIRS)}, and those with more than {MIN_MATCHES} matches)")

    centers = {im.name: im.projection_center() for im in rec.images.values()}
    rotations = {im.name: im.cam_from_world().rotation.matrix() for im in rec.images.values()}
    names = sorted(centers)
    scene = dict(camera_model=camera.cam.model.name, camera_params=np.asarray(camera.cam.params),
                 image_size=np.array([camera.cam.width, camera.cam.height]), principal_point=camera.center,
                 K=camera.K, image_names=names, camera_centers=np.array([centers[k] for k in names]),
                 camera_rotations=np.array([rotations[k] for k in names]), points3D=X, plane_inliers=inliers,
                 plane=np.r_[n, d], reprojection_error=float(all_errors.mean()))
    save_matches(Path(args.output), pairs, scene)
    print(f"\nWrote the matches and the ground truth of {len(pairs)} pairs to {args.output}")


if __name__ == "__main__":
    main()
