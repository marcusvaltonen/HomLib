"""
Figures 7 and 8 of Wadenbäck et al. (3DV 2026) on the Grossmünster Church, from the matches file of
grossmunster_ground_truth.py:
- Fig. 7: the 3D reconstruction, with the points of the plane of the wall (green), the other points
  (grey), and the camera positions of the images of the benchmark's pairs (with their optical axes);
- Fig. 8: (top) the inliers found by the two-sided equal solver (Wadenbäck et al., k1 = k2, in LOMSAC
  as in grossmunster_benchmark.py) on the pair DSC_2518 / DSC_2517; (bottom) the two images
  rectified with the estimated distortion: each pixel u of a rectified image, centered on the
  principal point, shows the pixel p = 2 u / (1 + sqrt(1 - 4 k |u|^2)) of the image, the inverse of
  the division model's undistortion u = p / (1 + k |p|^2). Lines that are straight in the scene
  should appear straight in them. The rectified images show the part of the image within its edges'
  midpoints, since the corners of a fisheye image go far off.

Requirements
------------
pip install matplotlib opencv-python

Usage
-----
    python grossmunster_figures.py data/grossmunster_sift.npz --data_dir data/fisheye_grossmunster --output_dir figures
"""

import argparse
import sys
from pathlib import Path

import cv2
import matplotlib
import numpy as np

import homlib
from grossmunster_data import FIGURE_PAIR, PAIRS, image_name, load_matches, pair_key

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import inlier_indices, make_estimator  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


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


def figure7(scene, numbers, path, elev=25, azim=-20):
    F = scene_frame(scene)
    X = scene["points3D"] @ F.T
    on_plane = scene["plane_inliers"].astype(bool)
    index = {name: i for i, name in enumerate(scene["image_names"])}
    ids = [index[image_name(n)] for n in numbers]
    C = scene["camera_centers"][ids] @ F.T
    axes = scene["camera_rotations"][ids][:, 2, :] @ F.T   # row 2 of R: the optical axis
    # Only the points around the wall and the cameras: the far ones would squeeze the scene
    lo = np.minimum(np.percentile(X[on_plane], 1, axis=0), C.min(0))
    hi = np.maximum(np.percentile(X[on_plane], 99, axis=0), C.max(0))
    margin = 0.1 * (hi - lo).max()
    near = np.all((X > lo - margin) & (X < hi + margin), axis=1)
    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(projection="3d")
    ax.scatter(*X[near & ~on_plane].T, s=0.2, c="0.65", depthshade=False)
    ax.scatter(*X[near & on_plane].T, s=0.4, c="tab:green", depthshade=False)
    ax.quiver(*C.T, *(axes * 0.1 * (hi - lo).max()).T, color="tab:red", linewidth=1)
    ax.scatter(*C.T, s=10, c="tab:red", depthshade=False)
    for n, c in zip(numbers, C):
        ax.text(*c, f" {n % 100}", fontsize=7)
    ax.set_box_aspect(np.ptp(np.vstack([X[near], C]), axis=0), zoom=1.3)
    ax.view_init(elev=elev, azim=azim)
    ax.set_axis_off()
    ax.set_title(f"The points of the plane (green) and the {len(numbers)} cameras (DSC_25xx)")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def rectify(image, center, k, width):
    """The image rectified with the division model's distortion k, in an image of this width: the
    part within the edges' midpoints of the image, at their undistorted positions"""
    h, w = image.shape[:2]
    # The undistorted edges' midpoints, around the principal point

    def undistort(p):
        return p / (1 + k * np.sum(p ** 2))
    left, right = undistort(np.array([-center[0], 0]))[0], undistort(np.array([w - 1 - center[0], 0]))[0]
    top, bottom = undistort(np.array([0, -center[1]]))[1], undistort(np.array([0, h - 1 - center[1]]))[1]
    if not (left < 0 < right and top < 0 < bottom):
        raise ValueError(f"the distortion k = {k} has no undistorted image of the edges' midpoints")
    height = int(round(width * (bottom - top) / (right - left)))
    u = np.stack(np.meshgrid(np.linspace(left, right, width), np.linspace(top, bottom, height)), -1)
    p = u * (2 / (1 + np.sqrt(np.maximum(1 - 4 * k * np.sum(u ** 2, -1), 0))))[..., None]
    maps = (p + center).astype(np.float32)
    # remap interpolates bilinearly: smooth the image first, which the output samples about 3 times
    # more coarsely
    blurred = cv2.GaussianBlur(image, (0, 0), 0.5 * w / width)
    return cv2.remap(blurred, maps[..., 0], maps[..., 1], cv2.INTER_LINEAR, borderMode=cv2.BORDER_CONSTANT)


def figure8(pairs, scene, image_dir, threshold, seed, path):
    a, b = FIGURE_PAIR
    key = pair_key(*sorted(FIGURE_PAIR))
    pair = pairs[key]
    x, y = (pair["x"], pair["y"]) if (a, b) == tuple(sorted(FIGURE_PAIR)) else (pair["y"], pair["x"])
    solver = homlib.Solver(homlib.SolverId.Wadenback3DV2026TwoSidedEqual)
    estimator = make_estimator("lomsac", solver, threshold, None, seed, min_iterations=500, final_least_squares=True)
    pose, stats = estimator.estimate(np.ascontiguousarray(x), np.ascontiguousarray(y))
    inliers = np.asarray(inlier_indices(stats))
    k = pose.distortion_parameter
    f = scene["K"][0, 0]
    print(f"Fig. 8: {len(inliers)} inliers of {x.shape[1]} matches, k = {k:.4g} (k f^2 = {k * f * f:.3f})")
    center = scene["principal_point"]
    images = [cv2.cvtColor(cv2.imread(str(image_dir / image_name(n))), cv2.COLOR_BGR2RGB) for n in (a, b)]
    w = images[0].shape[1]
    fig, axs = plt.subplot_mosaic([["top", "top"], ["left", "right"]], figsize=(12, 8.5))
    top = axs["top"]
    top.imshow(np.hstack(images))
    p1, p2 = x[:, inliers].T + center, y[:, inliers].T + center + [w, 0]
    for q1, q2 in zip(p1, p2):
        top.plot([q1[0], q2[0]], [q1[1], q2[1]], "-", color="lime", linewidth=0.8)
    top.plot(*p1.T, ".", color="yellow", markersize=3)
    top.plot(*p2.T, ".", color="yellow", markersize=3)
    top.set_title(f"DSC_{a} / DSC_{b}: the {len(inliers)} inliers of the two-sided equal solver")
    top.set_axis_off()
    for ax, n, image in zip((axs["left"], axs["right"]), (a, b), images):
        ax.imshow(rectify(image, center, k, 1600))
        ax.set_title(f"DSC_{n} rectified")
        ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Figures 7 and 8 of Wadenbäck et al. on the Grossmünster Church")
    parser.add_argument("matches", help="The matches file (grossmunster_ground_truth.py)")
    parser.add_argument("--data_dir", required=True, help="RadialSfM's fisheye_grossmunster (with images/)")
    parser.add_argument("--output_dir", default=".", help="Directory of the figures (default: .)")
    parser.add_argument("--threshold", type=float, default=3.0, help="Inlier threshold in pixels (default: 3)")
    parser.add_argument("--seed", type=int, default=1, help="Random seed of LOMSAC (default: 1)")
    args = parser.parse_args()
    pairs, scene = load_matches(Path(args.matches))
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    numbers = sorted({n for pair in PAIRS for n in pair})
    figure7(scene, numbers, out / "grossmunster_fig7.png")
    figure8(pairs, scene, Path(args.data_dir) / "images", args.threshold, args.seed, out / "grossmunster_fig8.png")
    print(f"Wrote {out / 'grossmunster_fig7.png'} and {out / 'grossmunster_fig8.png'}")


if __name__ == "__main__":
    main()
