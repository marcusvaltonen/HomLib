from __future__ import annotations

from homlib._core import (
    __doc__,
    PoseData,
    estimate_fitzgibbon_cvpr_2001_one_sided,
    estimate_fitzgibbon_cvpr_2001_two_sided_equal,
    estimate_kukelova_cvpr_2015_two_sided,
    estimate_kukelova_cvpr_2015_two_sided_6pt,
    estimate_nakano_icpr_2025_one_sided,
    estimate_wadenback_2025_one_sided,
    estimate_wadenback_2025_two_sided_equal,
    estimate_wadenback_2025_two_sided,
    lomsac_fitzgibbon_cvpr_2001_one_sided,
    lomsac_fitzgibbon_cvpr_2001_two_sided_equal,
    lomsac_kukelova_cvpr_2015_two_sided,
    lomsac_kukelova_cvpr_2015_two_sided_equal,
    lomsac_kukelova_cvpr_2015_two_sided_equal_6pt,
    lomsac_kukelova_cvpr_2015_two_sided_6pt,
    lomsac_nakano_icpr_2025_one_sided,
    lomsac_wadenback_2025_unsided,
    lomsac_wadenback_2025_one_sided,
    lomsac_wadenback_2025_two_sided_equal,
    lomsac_wadenback_2025_two_sided,
    LORansacOptions,
    RansacStatistics,
)
from homlib.version import __version__

__all__ = [
    "__doc__",
    "__version__",
    "PoseData",
    "estimate_fitzgibbon_cvpr_2001_one_sided",
    "estimate_fitzgibbon_cvpr_2001_two_sided_equal",
    "estimate_kukelova_cvpr_2015_two_sided",
    "estimate_kukelova_cvpr_2015_two_sided_6pt",
    "estimate_nakano_icpr_2025_one_sided",
    "estimate_wadenback_2025_one_sided",
    "estimate_wadenback_2025_two_sided_equal",
    "estimate_wadenback_2025_two_sided",
    "lomsac_fitzgibbon_cvpr_2001_one_sided",
    "lomsac_fitzgibbon_cvpr_2001_two_sided_equal",
    "lomsac_kukelova_cvpr_2015_two_sided",
    "lomsac_kukelova_cvpr_2015_two_sided_equal",
    "lomsac_kukelova_cvpr_2015_two_sided_equal_6pt",
    "lomsac_kukelova_cvpr_2015_two_sided_6pt",
    "lomsac_nakano_icpr_2025_one_sided",
    "lomsac_wadenback_2025_unsided",
    "lomsac_wadenback_2025_one_sided",
    "lomsac_wadenback_2025_two_sided_equal",
    "lomsac_wadenback_2025_two_sided",
    "LORansacOptions",
    "RansacStatistics",
]

import sys
import types
import time

from gluefactory.robust_estimators.base_estimator import BaseEstimator


class HomlibHomographyEstimator(BaseEstimator):
    # copied from poselib estimator
    default_conf = {"ransac_th": 2.0, "options": {}}

    def __init__(self, *args, **kwargs):
        self._image_size0 = None
        self._image_size1 = None
        self.solver = None
        super().__init__(*args, **kwargs)
        self.estimator = None

    def _init(self, conf):
        self._image_size0 = conf["image_size0"]
        self._image_size1 = conf["image_size1"]
        self.solver = getattr(sys.modules[__name__], conf["solver"])
        self.num_ransac_iter = conf["num_ransac_iter"]

    def _forward(self, data):
        if self._image_size0 is None or self._image_size1 is None:
            raise RuntimeError('Image size unknown. Cannot correct for distortion center assumed to be at the origin.')

        import torch
        import numpy as np
        pts0, pts1 = data["m_kpts0"], data["m_kpts1"]
        N = pts0.shape[0]

        options = LORansacOptions()
        options.squared_inlier_threshold = self.conf.ransac_th ** 2
        options.final_least_squares = True
        options.min_num_iterations = self.num_ransac_iter
        options.max_num_iterations = self.num_ransac_iter
        options.lo_starting_iterations = 20 # int(self.num_ransac_iter / 4)
        options.num_lsq_iterations = int(self.num_ransac_iter / 10)

        # Shift coordinates such that distortion parameters are located at the origin
        h1, w1 = self._image_size0
        h2, w2 = self._image_size1
        src_pts = pts0.numpy().mT - np.array([[h1, w1]]).T / 2.0
        dst_pts = pts1.numpy().mT - np.array([[h2, w2]]).T / 2.0
        T1 = np.eye(3)
        T1[0, 2] = -h1 / 2.0
        T1[1, 2] = -w1 / 2.0
        T2inv = np.eye(3)
        T2inv[0, 2] = h2 / 2.0
        T2inv[1, 2] = w2 / 2.0

        # Distort synthetically
        # Distort points synthetically (according to model)
        # kappa_gt = -1.0e-7
        # src_pts_dist = self._distort_points(src_pts, kappa_gt)
        # dst_pts_dist = self._distort_points(dst_pts, kappa_gt)

        pose_data, stats = self.solver(
            src_pts,
            dst_pts,
            options
        )

        """
        pose_data2, stats2 = lomsac_wadenback_2025_unsided(
            src_pts,
            dst_pts,
            options
        )

        print(f'one-sided: {stats}')
        print(f'unsided: {stats2}')
        print(f'diff norm {np.linalg.norm(pose_data.homography-pose_data2.homography)}')
        """
        if stats.best_num_inliers == 0:
            print("homlib failed to estimate homography")
            M = torch.eye(3, device=pts0.device, dtype=pts0.dtype)
            inl = torch.zeros_like(pts0[:, 0]).bool()
            success = False
        else:
            M = T2inv @ pose_data.homography @ T1
            M = torch.tensor(M).to(pts0)
            inl = torch.zeros(N, dtype=torch.bool)
            inl[torch.from_numpy(np.array(stats.inlier_indices))] = True
            success = True

        return {
            "success": success,
            "M_0to1": M,
            "inliers": inl,
            "time": stats.time,
        }

    def _distort_points(self, src_pts, kappa):
        """Use OpenCV to distort points using the one parameter division model."""
        import cv2
        import numpy as np
        x = np.vstack((src_pts, np.ones((1, src_pts.shape[1]), dtype=np.float32)))
        distortion_coeffs = np.array([0.0, 0.0, 0.0, 0.0, 0.0, kappa, 0.0, 0.0], dtype=np.float32)
        proj_d, _ = cv2.projectPoints(x.T, np.eye(3), np.zeros(3), np.eye(3), distortion_coeffs)
        proj_d = proj_d.T.squeeze()
        return proj_d


def register_homlib_estimator(
        target='gluefactory.robust_estimators.homography.homlib',
        source_type=HomlibHomographyEstimator
):
    # Create a new module with the GlueFactory-expected name
    proxy = types.ModuleType(target)

    # Rebind the class' __module__ so GlueFactory's loader filter matches
    source_type.__module__ = target
    setattr(proxy, source_type.__name__, source_type)

    # Install the proxy into sys.modules
    sys.modules[target] = proxy


# Call once before GlueFactory attempts to load the estimator
register_homlib_estimator()
