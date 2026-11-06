import unittest

import numpy as np

import homlib


def rot(rng=None, yaw_deg=None, angle_deg=10.0):
    """A rotation by yaw_deg about the vertical axis, or by angle_deg about a random axis"""
    if yaw_deg is not None:
        a = np.deg2rad(yaw_deg)
        return np.array([[np.cos(a), 0, np.sin(a)], [0, 1, 0], [-np.sin(a), 0, np.cos(a)]])
    ax = rng.normal(size=3)
    ax /= np.linalg.norm(ax)
    a = np.deg2rad(angle_deg)
    K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
    return np.eye(3) + np.sin(a) * K + (1 - np.cos(a)) * K @ K


def distort(u, k):
    """The division model's distorted points (k for these coordinates)"""
    return u * (2 / (1 + np.sqrt(1 - 4 * k * np.sum(u**2, 0))))


def matches(y, k_y, f_y, R, k_x, f_x):
    """The points x of image 1 matching the points y of image 2 (in pixels from their centers), for
    two cameras at the same viewpoint, with R from image 2's frame to image 1's; k for the
    coordinates divided by each image's focal length"""
    yn = y / f_y
    P = R @ np.vstack([yn / (1 + k_y * np.sum(yn**2, 0)), np.ones(y.shape[1])])
    return f_x * distort(P[:2] / P[2], k_x)


def minimal_samples(solver, x, y, scales, k1, k2, m=5, n=200):
    """Over n minimal samples of noise-free matches: the mean number of solutions, and the share of
    samples whose true solution (k1, k2 for the coordinates divided by the focal lengths) the
    solver keeps, the points of each sample divided by scales(idx) = (s1, s2, f1, f2)"""
    solver.normalize_image_coord = False   # the normalization is done here
    rng = np.random.default_rng(3)
    n_sols = n_kept = 0
    for _ in range(n):
        idx = rng.choice(x.shape[1], m, replace=False)
        s1, s2, f1, f2 = scales(idx)
        a, b = k1 * (s1 / f1)**2, k2 * (s2 / f2)**2   # the true k in these coordinates
        sols = solver.estimate(x[:, idx] / s1, y[:, idx] / s2)
        n_sols += len(sols)
        n_kept += any(abs(p.distortion_parameter - a) <= 1e-3 * abs(a) and
                      abs(p.distortion_parameter2 - b) <= 1e-3 * abs(b) for p in sols)
    return n_sols / n, n_kept / n


def f0(p):
    """The scale of the robust estimators' normalization: mean(|p|)/sqrt(2)"""
    return np.mean(np.linalg.norm(p, axis=0)) / np.sqrt(2)


class NormalizationTestCase(unittest.TestCase):

    def test_panorama(self):
        # A very large panorama (focal length 16000 px, slightly distorted) and a 1000 x 1000 view of
        # it from the same viewpoint (focal length 1000 px, k = -0.3). The estimators' normalization
        # takes its scale from the first image's points, here the panorama's, about 16 times the
        # view's: the view's true k is outside [-1, 0.1] in (almost) every sample. With each image divided by
        # its focal length (and normalize_image_coord off), it is -0.3, well inside.
        f_pano, f_view, k_pano, k_view = 16000.0, 1000.0, -0.1, -0.3
        y = np.random.default_rng(2).uniform(-500, 500, (2, 300))
        x = matches(y, k_view, f_view, rot(yaw_deg=2.0), k_pano, f_pano)
        for solver_class in (homlib.Wadenback3DV2026TwoSided, homlib.KukelovaCVPR2015TwoSided):
            solver = solver_class()
            solver.set_k_range(-1.0, 0.1)
            _, kept = minimal_samples(solver, x, y, lambda idx: (f0(x[:, idx]),) * 2 + (f_pano, f_view), k_pano, k_view)
            self.assertLess(kept, 0.05)
            _, kept = minimal_samples(solver, x, y, lambda idx: (f_pano, f_view, f_pano, f_view), k_pano, k_view)
            self.assertGreater(kept, 0.95)

    def test_tight_range_with_known_focal_length(self):
        # The same camera (4000 x 3000 px, a 14 mm lens in 35 mm equivalent, k = -0.3 for the
        # coordinates divided by its focal length) in two views rotated by 10 degrees. With the
        # points divided by the focal length, a tight range of real lenses, [-0.6, 0.05], keeps the
        # true solution and leaves fewer solutions per sample than the default [-1, 0.1] in the
        # estimators' normalized coordinates (about [-3, 0.3] in the focal-normalized ones here).
        f, k = 14 / 43.27 * 5000, -0.3
        rng = np.random.default_rng(1)
        x = np.c_[rng.uniform(-2000, 2000, 3000), rng.uniform(-1500, 1500, 3000)].T
        y = matches(x, k, f, rot(rng), k, f)
        inside = (np.abs(y[0]) < 2000) & (np.abs(y[1]) < 1500)
        x, y = x[:, inside][:, :500], y[:, inside][:, :500]
        solver = homlib.Wadenback3DV2026TwoSided()
        solver.set_k_range(-1.0, 0.1)
        n_default, kept_default = minimal_samples(solver, x, y, lambda idx: (f0(x[:, idx]),) * 2 + (f, f), k, k)
        solver.set_k_range(-0.6, 0.05)
        n_tight, kept_tight = minimal_samples(solver, x, y, lambda idx: (f, f, f, f), k, k)
        self.assertGreater(kept_default, 0.95)
        self.assertGreater(kept_tight, 0.95)
        self.assertLess(n_tight, 0.6 * n_default)


if __name__ == '__main__':
    unittest.main()
