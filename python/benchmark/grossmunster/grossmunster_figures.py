"""
Figures 7 and 8 of Wadenbäck et al. (3DV 2026) on the Grossmünster Church, from the matches file of
grossmunster_ground_truth.py:
- Fig. 7: the 3D reconstruction seen from in front of the wall, as COLMAP's viewer shows it, with
  the points of the plane of the wall (green), the other points (black), a grid on the ground, and
  the cameras of the benchmark's images as frustums showing their images;
- Fig. 8: (top) the inliers found by the two-sided equal solver (Wadenbäck et al., k1 = k2, in LOMSAC
  as in grossmunster_benchmark.py) on the pair DSC_2518 / DSC_2517, as yellow lines; (bottom) the
  images DSC_2512, DSC_2514 and DSC_2519, as in the paper, each rectified with the distortion k that
  the solver estimates on its pair (RECTIFIED_PAIRS in grossmunster_data.py): each pixel u of a
  rectified image, centered on the principal point, shows the pixel p = 2 u / (1 + sqrt(1 - 4 k |u|^2))
  of the image, the inverse of the division model's undistortion u = p / (1 + k |p|^2). Lines that are
  straight in the scene should appear straight in them. The whole image is rectified: its corners go
  far off, which gives the star shape.

Requirements
------------
pip install opencv-python

Usage
-----
    python grossmunster_figures.py data/grossmunster_sift.npz --data_dir data/fisheye_grossmunster --output_dir figures
"""

import argparse
import sys
from pathlib import Path

import cv2
import numpy as np

import homlib
from grossmunster_data import (FIGURE_PAIR, LOMSAC_SETTINGS, MIN_ITERATIONS, PAIRS, RECTIFIED_PAIRS, THRESHOLD,
                               image_name, load_matches, pair_key)

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import inlier_indices, make_estimator  # noqa: E402


def scene_frame(scene):
    """The rows of a rotation from the world to a frame of the wall: along the wall (the direction
    in which the cameras move), its normal (towards the cameras), and up (towards most of the points
    of the plane, above the cameras)"""
    n = scene["plane"][:3]
    C = scene["camera_centers"]
    if np.mean(C @ n - scene["plane"][3]) < 0:
        n = -n
    along = np.linalg.svd(C - C.mean(0))[2][0]
    along -= (along @ n) * n
    along /= np.linalg.norm(along)
    up = np.cross(along, n)   # a right-handed frame: along x n = up
    if np.median(scene["points3D"][scene["plane_inliers"].astype(bool)] @ up) < np.median(C @ up):
        along, up = -along, -up
    return np.array([along, n, up])


def look_at(eye, target, up):
    """The rotation (rows: right, down, forward) and center of a virtual camera at eye looking at target"""
    forward = (target - eye) / np.linalg.norm(target - eye)
    right = np.cross(forward, up)
    right /= np.linalg.norm(right)
    return np.array([right, np.cross(forward, right), forward]), eye


def figure7(scene, numbers, image_dir, path, size=(1100, 800), thumb_width=200):
    """The reconstruction seen from in front of the wall, as COLMAP's viewer shows it: the points of
    the plane in green, the others in black, a grid on the ground, and the cameras of the benchmark's
    images as frustums with their images"""
    F = scene_frame(scene)            # rows: along the wall, its normal (towards the cameras), up
    X = scene["points3D"] @ F.T
    on_plane = scene["plane_inliers"].astype(bool)
    index = {name: i for i, name in enumerate(scene["image_names"])}
    ids = [index[image_name(n)] for n in numbers]
    C = scene["camera_centers"][ids] @ F.T
    R = scene["camera_rotations"][ids] @ F.T   # rows: the cameras' axes in the wall's frame
    wall = X[on_plane]
    lo, hi = np.percentile(wall, 2, axis=0), np.percentile(wall, 98, axis=0)
    extent = (hi - lo).max()
    ground = np.percentile(X[:, 2], 1)
    target = np.array([(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2])
    eye = target + np.array([0.35, 2.1, 0.55]) * extent
    Rv, cv = look_at(eye, target, np.array([0, 0, 1.0]))
    w, h = size
    f = 0.9 * w

    def project(P):
        Q = (np.atleast_2d(P) - cv) @ Rv.T
        return Q[:, :2] / Q[:, 2:] * f + [w / 2, h / 2], Q[:, 2]

    canvas = np.full((h, w, 3), 255, np.uint8)
    # The grid on the ground
    step = extent / 8
    for v in np.arange(-6, 7) * step:
        for a, b in (([lo[0] + v * 1.5, -4 * extent, ground], [lo[0] + v * 1.5, 4 * extent, ground]),
                     ([-4 * extent, target[1] + v * 1.5, ground], [4 * extent, target[1] + v * 1.5, ground])):
            t = np.linspace(0, 1, 200)[:, None]
            q, z = project(np.array(a) + t * (np.array(b) - np.array(a)))
            q = q[z > 0.05 * extent]
            cv2.polylines(canvas, [np.round(q * 16).astype(np.int32)], False, (225, 225, 225), 1, cv2.LINE_AA, 4)
    # The points, far ones first
    q, z = project(X)
    for keep, color in ((~on_plane, (40, 40, 40)), (on_plane, (30, 200, 30))):
        for (u, v) in np.round(q[keep & (z > 0)]).astype(int):
            if 0 <= u < w and 0 <= v < h:
                cv2.circle(canvas, (u, v), 1, color, -1, cv2.LINE_AA)
    # The cameras, far ones first, as frustums whose base shows the image
    depth = 0.1 * extent
    order = np.argsort(-project(C)[1])
    for i in order:
        n = numbers[i]
        image = cv2.imread(str(image_dir / image_name(n)))
        image = cv2.resize(image, (thumb_width, int(thumb_width * image.shape[0] / image.shape[1])))
        ih, iw = image.shape[:2]
        corners = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]]) * [0.75, 0.5]
        base = C[i] + (np.c_[corners, np.ones(4)] * depth) @ R[i]
        qb, zb = project(base)
        qc, zc = project(C[i])
        if np.any(zb <= 0) or zc[0] <= 0:
            continue
        H = cv2.getPerspectiveTransform(np.float32([[0, 0], [iw, 0], [iw, ih], [0, ih]]), np.float32(qb))
        warped = cv2.warpPerspective(image[:, :, ::-1], H, (w, h))
        mask = cv2.warpPerspective(np.full((ih, iw), 255, np.uint8), H, (w, h)) > 0
        canvas[mask] = (0.85 * warped[mask] + 0.15 * 255).astype(np.uint8)
        frame = np.round(qb * 16).astype(np.int32)
        apex = np.round(qc[0] * 16).astype(np.int32)
        cv2.polylines(canvas, [frame], True, (120, 140, 150), 1, cv2.LINE_AA, 4)
        for c in frame:
            cv2.line(canvas, tuple(apex), tuple(c), (120, 140, 150), 1, cv2.LINE_AA, 4)
    cv2.imwrite(str(path), canvas[:, :, ::-1])


def rectify(image, center, k, width):
    """The whole image rectified with the division model's distortion k, in an image of this width
    whose corners are the undistorted corners of the image (or 4 times the undistorted edges'
    midpoints, if the corners go farther), on a white background"""
    h, w = image.shape[:2]

    def undistort(p):
        return p / (1 + k * np.sum(p ** 2, -1, keepdims=True))
    corners = np.array([[0, 0], [w - 1, 0], [w - 1, h - 1], [0, h - 1]]) - center
    midpoints = np.array([[0, center[1]], [w - 1, center[1]], [center[0], 0], [center[0], h - 1]]) - center
    limit = 4 * np.abs(undistort(midpoints)).max(0)
    denominator = 1 + k * np.sum(corners ** 2, -1)
    half = np.minimum(np.abs(undistort(corners)).max(0), limit) if np.all(denominator > 0) else limit
    height = int(round(width * half[1] / half[0]))
    u = np.stack(np.meshgrid(np.linspace(-half[0], half[0], width), np.linspace(-half[1], half[1], height)), -1)
    p = u * (2 / (1 + np.sqrt(np.maximum(1 - 4 * k * np.sum(u ** 2, -1), 0))))[..., None]
    maps = (p + center).astype(np.float32)
    # remap interpolates bilinearly: smooth the image first (lightly: the output samples its centre about
    # 10 times more coarsely, and its edges more finely)
    blurred = cv2.GaussianBlur(image, (0, 0), max(0.25 * 2 * half[0] / width, 0.5))
    return cv2.remap(blurred, maps[..., 0], maps[..., 1], cv2.INTER_LINEAR, borderMode=cv2.BORDER_CONSTANT,
                     borderValue=(255, 255, 255))


def estimate(pairs, a, b, threshold, seed):
    """The two-sided equal solver in LOMSAC on the pair (a, b): (pose, inliers, x, y), with x in image a"""
    pair = pairs[pair_key(*sorted((a, b)))]
    x, y = (pair["x"], pair["y"]) if a < b else (pair["y"], pair["x"])
    solver = homlib.Solver(homlib.SolverId.Wadenback3DV2026TwoSidedEqual)
    estimator = make_estimator("lomsac", solver, threshold, None, seed,
                               min_iterations=MIN_ITERATIONS, **LOMSAC_SETTINGS)
    pose, stats = estimator.estimate(np.ascontiguousarray(x), np.ascontiguousarray(y))
    return pose, np.asarray(inlier_indices(stats)), x, y


def figure8(pairs, scene, image_dir, threshold, seed, path, height=680):
    """(Top) the inliers of the two-sided equal solver on the pair of Fig. 8, as yellow lines; (bottom)
    three images rectified with the distortion it estimates on their pair (RECTIFIED_PAIRS)"""
    a, b = FIGURE_PAIR
    pose, inliers, x, y = estimate(pairs, a, b, threshold, seed)
    center = scene["principal_point"]
    print(f"Fig. 8: {len(inliers)} inliers of {x.shape[1]} matches, k = {pose.distortion_parameter:.4g}")
    images = [cv2.imread(str(image_dir / image_name(n))) for n in (a, b)]
    scale = height / images[0].shape[0]
    top = np.hstack([cv2.resize(im, None, fx=scale, fy=scale, interpolation=cv2.INTER_AREA) for im in images])
    w = top.shape[1] // 2
    p1, p2 = (x[:, inliers].T + center) * scale, (y[:, inliers].T + center) * scale + [w, 0]
    lines = top.copy()
    for q1, q2 in zip(p1, p2):
        cv2.line(lines, tuple(np.round(q1 * 16).astype(int)), tuple(np.round(q2 * 16).astype(int)), (0, 255, 255), 1,
                 cv2.LINE_AA, 4)
    top = cv2.addWeighted(lines, 0.6, top, 0.4, 0)   # semi-transparent lines, as the paper's
    rectified = []
    for n, m in RECTIFIED_PAIRS:
        pose_n = estimate(pairs, n, m, threshold, seed)[0]
        image = cv2.imread(str(image_dir / image_name(n)))
        rectified.append(rectify(image, center, pose_n.distortion_parameter, top.shape[1] // 3 - 20))
        print(f"Fig. 8: DSC_{n} rectified with k = {pose_n.distortion_parameter:.4g} (pair {n}-{m})")
    rh = max(r.shape[0] for r in rectified)
    row = np.full((rh, top.shape[1], 3), 255, np.uint8)
    for i, r in enumerate(rectified):
        x0 = i * (top.shape[1] // 3) + 10
        row[(rh - r.shape[0]) // 2:(rh - r.shape[0]) // 2 + r.shape[0], x0:x0 + r.shape[1]] = r
    gap = np.full((20, top.shape[1], 3), 255, np.uint8)
    cv2.imwrite(str(path), np.vstack([top, gap, row]))


def main():
    parser = argparse.ArgumentParser(description="Figures 7 and 8 of Wadenbäck et al. on the Grossmünster Church")
    parser.add_argument("matches", help="The matches file (grossmunster_ground_truth.py)")
    parser.add_argument("--data_dir", required=True, help="RadialSfM's fisheye_grossmunster (with images/)")
    parser.add_argument("--output_dir", default=".", help="Directory of the figures (default: .)")
    parser.add_argument("--threshold", type=float, default=THRESHOLD,
                        help="Inlier threshold in pixels (default: %(default)s)")
    parser.add_argument("--seed", type=int, default=1, help="Random seed of LOMSAC (default: 1)")
    args = parser.parse_args()
    pairs, scene = load_matches(Path(args.matches))
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    numbers = sorted({n for pair in PAIRS for n in pair})
    figure7(scene, numbers, Path(args.data_dir) / "images", out / "grossmunster_fig7.png")
    figure8(pairs, scene, Path(args.data_dir) / "images", args.threshold, args.seed, out / "grossmunster_fig8.png")
    print(f"Wrote {out / 'grossmunster_fig7.png'} and {out / 'grossmunster_fig8.png'}")


if __name__ == "__main__":
    main()
