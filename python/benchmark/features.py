"""Feature type classes for the benchmarks: extraction + matching, named by feature type.

Each class encapsulates a feature detector, descriptor, and matcher.  The base class
FeatureMatcher (in common.py) provides a concrete match() that calls extract() on each image
and _match_features() on the results.  End-to-end systems (LoMa, ALoMa) override match().

Feature types with a single natural matcher implement _match_features() directly.  Types with
multiple compatible matchers (SuperPoint: LightGlue or SuperGlue) accept a ``matcher``
constructor parameter and dispatch in _match_features().
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

import homlib
from common import (
    AffineMatchResult,
    FeatureMatcher,
    FeatureRegistry,
    FeatureUpgrader,
    OrientationMatchResult,
    PointMatchResult,
    ScaleOrientationMatchResult,
)
from matchers import BFMatcher, LightGlueMatcher


# ---------------------------------------------------------------------------
# SIFT
# ---------------------------------------------------------------------------

class SIFT(FeatureMatcher):
    """OpenCV SIFT, with BFMatcher + Lowe's ratio test (default) or LightGlue."""
    supported_features = {homlib.FeatureType.Points, homlib.FeatureType.Orientation}
    cli_matchers = {"bfmatcher", "lightglue"}

    @classmethod
    def from_args(cls, args):
        if args.matcher == "lightglue":
            return cls(num_features=args.num_features or 2048,
                       matcher="lightglue", resize=args.resize)
        return cls(num_features=args.num_features or 4000,
                   lowe_ratio=0.8 if args.lowe_ratio is None else args.lowe_ratio)

    def __init__(self, num_features=4000, lowe_ratio=0.8, matcher="bfmatcher", resize=None):
        import cv2
        super().__init__()
        self._matcher_name = matcher
        self._resize = resize
        if matcher == "bfmatcher":
            try:
                self._sift = cv2.SIFT_create(num_features, enable_precise_upscale=True)
            except (TypeError, cv2.error):
                self._sift = cv2.SIFT_create(num_features)
            self._matcher = BFMatcher(lowe_ratio=lowe_ratio)
        elif matcher == "lightglue":
            import torch
            from lightglue import SIFT as LightGlueSIFT
            self._device = ("cuda" if torch.cuda.is_available()
                            else "mps" if torch.backends.mps.is_available() else "cpu")
            self._extractor = LightGlueSIFT(
                max_num_keypoints=num_features, origin="center").eval().to(self._device)
            if self._extractor.conf.backend == "opencv":
                c = self._extractor.conf
                try:
                    self._extractor.sift = cv2.SIFT_create(
                        contrastThreshold=c.detection_threshold, nfeatures=c.max_num_keypoints,
                        edgeThreshold=c.edge_threshold, nOctaveLayers=c.num_octaves,
                        enable_precise_upscale=True)
                except (TypeError, cv2.error):
                    pass
            self._matcher = LightGlueMatcher("sift", self._device)
        else:
            raise ValueError(f"Unknown matcher: {matcher}")

    def _extract(self, path):
        if self._matcher_name == "bfmatcher":
            return self._extract_opencv(path)
        return self._extract_lightglue(path)

    def _extract_opencv(self, path):
        import cv2
        img = cv2.imread(str(path), cv2.IMREAD_GRAYSCALE)
        kps, descs = self._sift.detectAndCompute(img, None)
        return dict(keypoints=kps, descriptors=descs)

    def _extract_lightglue(self, path):
        from lightglue.utils import load_image
        options = dict(resize=self._resize, side="short") if self._resize else dict(resize=None)
        return self._extractor.extract(load_image(str(path)).to(self._device), **options)

    def _match_features(self, f1, f2):
        idx1, idx2, scores = self._matcher.match(f1, f2)
        if self._matcher_name == "bfmatcher":
            k1, k2 = f1["keypoints"], f2["keypoints"]
            return ScaleOrientationMatchResult(
                points1=np.array([k1[i].pt for i in idx1]).reshape(-1, 2),
                points2=np.array([k2[i].pt for i in idx2]).reshape(-1, 2),
                scores=scores,
                angles1=np.deg2rad(np.array([k1[i].angle for i in idx1])),
                angles2=np.deg2rad(np.array([k2[i].angle for i in idx2])),
                scales1=np.array([k1[i].size for i in idx1]),
                scales2=np.array([k2[i].size for i in idx2]),
            )
        from lightglue.utils import rbd
        f1, f2 = rbd(f1), rbd(f2)
        return PointMatchResult(
            points1=f1["keypoints"][idx1].cpu().numpy().astype(float),
            points2=f2["keypoints"][idx2].cpu().numpy().astype(float),
            scores=scores,
        )


# ---------------------------------------------------------------------------
# ORB
# ---------------------------------------------------------------------------

class ORB(FeatureMatcher):
    """OpenCV ORB with BFMatcher (Hamming) + Lowe's ratio test.  Orientation only, no meaningful scale."""
    supported_features = {homlib.FeatureType.Points, homlib.FeatureType.Orientation}

    def __init__(self, num_features=4000, lowe_ratio=0.8):
        import cv2
        super().__init__()
        self._orb = cv2.ORB_create(num_features)
        self._matcher = BFMatcher(lowe_ratio=lowe_ratio, norm=cv2.NORM_HAMMING)

    def _extract(self, path):
        import cv2
        img = cv2.imread(str(path), cv2.IMREAD_GRAYSCALE)
        kps, descs = self._orb.detectAndCompute(img, None)
        return dict(keypoints=kps, descriptors=descs)

    def _match_features(self, f1, f2):
        idx1, idx2, scores = self._matcher.match(f1, f2)
        k1, k2 = f1["keypoints"], f2["keypoints"]
        return OrientationMatchResult(
            points1=np.array([k1[i].pt for i in idx1]).reshape(-1, 2),
            points2=np.array([k2[i].pt for i in idx2]).reshape(-1, 2),
            scores=scores,
            angles1=np.deg2rad(np.array([k1[i].angle for i in idx1])),
            angles2=np.deg2rad(np.array([k2[i].angle for i in idx2])),
        )


# ---------------------------------------------------------------------------
# LightGlue-based features (SuperPoint, DISK, ALIKED, SIFT)
# ---------------------------------------------------------------------------

class _LightGlueFeature(FeatureMatcher):
    """Base for feature types matched with LightGlue.  The LightGlue extractor and matcher
    are selected by the lowercase class name (SuperPoint → "superpoint", etc.)."""
    supported_features = {homlib.FeatureType.Points}
    cli_matchers = {"lightglue"}

    @classmethod
    def from_args(cls, args):
        return cls(num_features=args.num_features or 2048, resize=args.resize)

    def __init__(self, num_features=2048, resize=None, matcher="lightglue"):
        import torch
        super().__init__()
        self._device = ("cuda" if torch.cuda.is_available()
                        else "mps" if torch.backends.mps.is_available() else "cpu")
        self._resize = resize
        self._matcher_name = matcher
        self._setup_extractor(num_features)
        self._setup_matcher()

    def _setup_extractor(self, num_features):
        from lightglue import ALIKED, DISK, SuperPoint
        cls = {"superpoint": SuperPoint, "disk": DISK, "aliked": ALIKED}[self.name]
        self._extractor = cls(max_num_keypoints=num_features, origin="center").eval().to(self._device)

    def _setup_matcher(self):
        if self._matcher_name == "lightglue":
            self._matcher = LightGlueMatcher(self.name, self._device)
        else:
            raise ValueError(f"Unknown matcher: {self._matcher_name}")

    def _extract(self, path):
        from lightglue.utils import load_image
        options = dict(resize=self._resize, side="short") if self._resize else dict(resize=None)
        return self._extractor.extract(load_image(str(path)).to(self._device), **options)

    def _match_features(self, f1, f2):
        from lightglue.utils import rbd
        idx1, idx2, scores = self._matcher.match(f1, f2)
        f1, f2 = rbd(f1), rbd(f2)
        return PointMatchResult(
            points1=f1["keypoints"][idx1].cpu().numpy().astype(float),
            points2=f2["keypoints"][idx2].cpu().numpy().astype(float),
            scores=scores,
        )


class SuperPoint(_LightGlueFeature):
    """SuperPoint features, matched with LightGlue (default) or SuperGlue."""

    def __init__(self, num_features=2048, resize=None, matcher="lightglue"):
        super().__init__(num_features=num_features, resize=resize, matcher=matcher)

    def _setup_matcher(self):
        if self._matcher_name == "superglue":
            raise NotImplementedError("SuperGlue matcher not yet implemented")
        super()._setup_matcher()

    def _match_features(self, features1, features2):
        if self._matcher_name == "superglue":
            return self._match_superglue(features1, features2)
        return super()._match_features(features1, features2)

    def _match_superglue(self, features1, features2):
        raise NotImplementedError("SuperGlue matcher not yet implemented")


class DISK(_LightGlueFeature):
    """DISK features matched with LightGlue."""

    def __init__(self, num_features=2048, resize=None):
        super().__init__(num_features=num_features, resize=resize)


class ALIKED(_LightGlueFeature):
    """ALIKED features matched with LightGlue."""

    def __init__(self, num_features=2048, resize=None):
        super().__init__(num_features=num_features, resize=resize)


# ---------------------------------------------------------------------------
# LoMa
# ---------------------------------------------------------------------------

class LoMa(FeatureMatcher):
    """LoMa (Nordström et al., ECCV 2026): end-to-end detect + describe + match."""
    supported_features = {homlib.FeatureType.Points}
    cli_matchers = set()

    @classmethod
    def from_args(cls, args):
        return cls(num_keypoints=args.num_features or 2048)

    def __init__(self, num_keypoints=2048):
        import torch
        from loma import LoMa as LoMaModel, LoMaB
        super().__init__()
        self._model = LoMaModel(LoMaB(num_keypoints=num_keypoints))
        self._torch = torch

    def _extract(self, path):
        raise NotImplementedError("LoMa is end-to-end: use match() directly")

    def _match_features(self, features1, features2):
        raise NotImplementedError("LoMa is end-to-end: use match() directly")

    def match(self, path1, path2):
        from loma.loma import filter_matches, to_pixel_coords
        try:
            with self._torch.inference_mode():
                k1, d1, h1, w1 = self._model.detect_and_describe(str(path1), None)
                k2, d2, h2, w2 = self._model.detect_and_describe(str(path2), None)
                m0, _, scores0, _ = filter_matches(
                    self._model(k1, k2, d1, d2)["scores"],
                    self._model.cfg.filter_threshold,
                )
            valid = m0[0] > -1
            kpts1 = to_pixel_coords(k1[0][self._torch.where(valid)[0]], h1, w1).cpu().numpy()
            kpts2 = to_pixel_coords(k2[0][m0[0][valid]], h2, w2).cpu().numpy()
            scores = scores0[0][valid].cpu().numpy()
        except Exception as e:
            print(f"LoMa failed on {path1} and {path2}: {e}")
            kpts1 = kpts2 = np.zeros((0, 2))
            scores = np.zeros(0)
        # LoMa's (0, 0) is the top-left corner: shift to OpenCV's convention
        return PointMatchResult(
            points1=np.asarray(kpts1, dtype=float).reshape(-1, 2) - 0.5,
            points2=np.asarray(kpts2, dtype=float).reshape(-1, 2) - 0.5,
            scores=np.asarray(scores, dtype=float).reshape(-1),
        )


# ---------------------------------------------------------------------------
# AffNet (kornia)
# ---------------------------------------------------------------------------

class AffNet(FeatureMatcher):
    """AffNet + HardNet (via kornia): SIFT detector with AffNet affine estimation and
    HardNet descriptors, matched with BFMatcher + Lowe's ratio test."""
    supported_features = {
        homlib.FeatureType.Points, homlib.FeatureType.Orientation, homlib.FeatureType.Affine,
    }

    def __init__(self, num_features=4000, lowe_ratio=0.8):
        import torch
        import kornia.feature as KF
        from kornia_moons.feature import OpenCVDetectorWithAffNetKornia
        import cv2

        super().__init__()
        self._device = "cuda" if torch.cuda.is_available() else "cpu"
        detector = OpenCVDetectorWithAffNetKornia(cv2.SIFT_create(num_features), max_kpts=num_features)
        descriptor = KF.LAFDescriptor(KF.HardNet(True)).eval()
        self._feature = KF.LocalFeature(detector, descriptor)
        self._matcher = BFMatcher(lowe_ratio=lowe_ratio)

    def _extract(self, path):
        import cv2
        import kornia as K
        import torch
        img = cv2.imread(str(path))
        img_t = K.image_to_tensor(img, False).float() / 255.0
        img_t = img_t.to(self._device)
        with torch.no_grad():
            lafs, _, descs = self._feature(img_t)
        return dict(
            lafs=lafs.squeeze(0).cpu().numpy(),
            descriptors=descs.squeeze(0).cpu().numpy(),
        )

    def _match_features(self, f1, f2):
        idx1, idx2, scores = self._matcher.match(f1, f2)
        return AffineMatchResult.from_lafs(f1["lafs"][idx1], f2["lafs"][idx2], scores=scores)


# ---------------------------------------------------------------------------
# DenseAffine
# ---------------------------------------------------------------------------

class DenseAffine(FeatureMatcher):
    """DenseAffine (Sun et al., CVPR 2025): end-to-end dense affine pipeline.
    DKMv3 -> EstiNet scale+angle -> HardNet -> affine correspondences."""
    supported_features = {
        homlib.FeatureType.Points, homlib.FeatureType.Orientation, homlib.FeatureType.Affine,
    }

    def __init__(self, repo_path, weights_path=None, device=None):
        super().__init__()
        self._repo_path = repo_path
        self._weights_path = weights_path
        self._device = device

    def _extract(self, path):
        raise NotImplementedError("DenseAffine is end-to-end: use match() directly")

    def _match_features(self, features1, features2):
        raise NotImplementedError("DenseAffine is end-to-end: use match() directly")

    def match(self, path1, path2):
        raise NotImplementedError(
            "DenseAffine requires the DenseAffine repo (https://github.com/stilcrad/DenseAffine). "
            "Integration pending."
        )


# ---------------------------------------------------------------------------
# POP
# ---------------------------------------------------------------------------

class POP(FeatureMatcher):
    """POP (elvintanhust/POP-Interest-Point): learned unsupervised interest point
    detector/descriptor.  Points only (no orientation or scale)."""
    supported_features = {homlib.FeatureType.Points}

    def __init__(self, model_path=None, num_features=1000, device=None):
        super().__init__()
        self._model_path = model_path
        self._num_features = num_features
        self._device = device

    def _extract(self, path):
        raise NotImplementedError(
            "POP requires the POP repo (https://github.com/elvintanhust/POP-Interest-Point). "
            "Integration pending."
        )

    def _match_features(self, features1, features2):
        raise NotImplementedError("POP integration pending")


# ---------------------------------------------------------------------------
# Upgrader stubs
# ---------------------------------------------------------------------------

class AffineDecompositionUpgrader(FeatureUpgrader):
    """Decompose per-keypoint 2x2 affine frames into angles and scales.

    Extracts angle = atan2(a21, a11) and scale = sqrt(|det(A)|) from each
    individual affine frame, producing a ScaleOrientationMatchResult suitable
    for orientation or scale+orientation solvers."""
    input_type = AffineMatchResult
    output_type = ScaleOrientationMatchResult

    def upgrade(self, result, **kwargs):
        a1, a2 = result.affines1, result.affines2
        return ScaleOrientationMatchResult(
            points1=result.points1, points2=result.points2, scores=result.scores,
            angles1=np.arctan2(a1[:, 1, 0], a1[:, 0, 0]),
            angles2=np.arctan2(a2[:, 1, 0], a2[:, 0, 0]),
            scales1=np.sqrt(np.abs(np.linalg.det(a1))),
            scales2=np.sqrt(np.abs(np.linalg.det(a2))),
        )


class BarathAffineUpgrader(FeatureUpgrader):
    """Upgrade scale+orientation matches to affine correspondences using the epipolar constraint.

    Given known SIFT angles (alpha1, alpha2) and scales (q1, q2) and a fundamental matrix F,
    solve for the unknown affine parameters (q_u, q_v, w) in closed form per correspondence.

    Reference: Barath, 'On Making SIFT Features Affine Covariant' (IJCV 2023, Eqs. 11-18).
    A = R2 @ U @ R1^T  where U = [[q_u, w], [0, q_v]] and q_u = q2^2 / (q1^2 * q_v).
    """
    input_type = ScaleOrientationMatchResult
    output_type = AffineMatchResult

    def upgrade(self, result, *, F):
        n = len(result.points1)
        affines1 = np.zeros((n, 2, 2))
        affines2 = np.zeros((n, 2, 2))

        for i in range(n):
            a1, a2 = result.angles1[i], result.angles2[i]
            q1, q2 = result.scales1[i], result.scales2[i]
            u1, v1 = result.points1[i]
            u2, v2 = result.points2[i]
            c1, s1 = np.cos(a1), np.sin(a1)
            c2, s2 = np.cos(a2), np.sin(a2)

            # f values from F (row-major): F = [[f1,f2,f3],[f4,f5,f6],[f7,f8,f9]]
            f1, f2, f3, f4, f5, f6, f7, f8, f9 = F.ravel()

            # Constants from Eq. 13-17 of the paper
            B = u1 * f1 + v1 * f2 + f3
            C = u1 * f4 + v1 * f5 + f6
            D = -u2 * f1 - v2 * f4 - f7
            E = -u2 * f2 - v2 * f5 - f8

            G = B * c2 * c1 + C * s2 * c1
            H = B * s2 * s1 - C * c2 * s1
            I_val = -B * c2 * s1 - C * s2 * s1
            J = B * c2 * s1 + C * s2 * s1

            # Solve for q_v from the linear equation (Eq. 17 simplified):
            # The coefficient of q_v^2 is K - G*H/I_val which is zero by trig identities.
            # So: (G*D/I_val - E) * q_v + J * q2^2/q1^2 - G^2 * q2^2 / (I_val * q1^2) = 0
            q_ratio_sq = q2 ** 2 / (q1 ** 2 + 1e-30)

            if abs(I_val) < 1e-15:
                # Degenerate: epipolar line nearly parallel to feature axis
                affines1[i] = q1 * np.array([[c1, -s1], [s1, c1]])
                affines2[i] = q2 * np.array([[c2, -s2], [s2, c2]])
                continue

            coeff_qv = G * D / I_val - E
            const_term = J * q_ratio_sq - G ** 2 * q_ratio_sq / I_val

            if abs(coeff_qv) < 1e-15:
                affines1[i] = q1 * np.array([[c1, -s1], [s1, c1]])
                affines2[i] = q2 * np.array([[c2, -s2], [s2, c2]])
                continue

            q_v = -const_term / coeff_qv
            q_u = q_ratio_sq / (q_v + 1e-30)
            w = (D - G * q_u - H * q_v) / I_val

            # Reconstruct A = R2 @ [[q_u, w], [0, q_v]] @ R1^T
            R1 = np.array([[c1, -s1], [s1, c1]])
            R2 = np.array([[c2, -s2], [s2, c2]])
            U = np.array([[q_u, w], [0, q_v]])
            A = R2 @ U @ R1.T

            # Individual affine frames: LAF1 = q1 * R1, LAF2 = A @ LAF1
            affines1[i] = q1 * R1
            affines2[i] = A @ affines1[i]

        return AffineMatchResult(
            points1=result.points1, points2=result.points2, scores=result.scores,
            angles1=result.angles1, angles2=result.angles2,
            scales1=result.scales1, scales2=result.scales2,
            affines1=affines1, affines2=affines2,
        )


class S3EstiUpgrader(FeatureUpgrader):
    """Upgrade point matches to scale+orientation using S3Esti (EstiNet).

    Estimates per-keypoint scale and orientation from image patches using the EstiNet CNN
    (300 scale bins, 360 angle bins).  Requires the S3Esti repo
    (thirdparty/S3Esti or a path to it) and model weights (S3Esti_ep30.pth).

    Reference: S3Esti — EstiNet with multi-scale patch extraction.
    """
    input_type = PointMatchResult
    output_type = ScaleOrientationMatchResult

    def __init__(self, s3esti_path, checkpoint="S3Esti_ep30.pth", device=None,
                 patch_size=32, scale_ratios=(0.5, 1, 2), num_patches=4):
        import sys as _sys
        import torch
        self._device = device or ("cuda" if torch.cuda.is_available() else "cpu")
        self._patch_size = patch_size
        self._scale_ratios = list(scale_ratios)
        self._num_patches = num_patches

        s3esti_path = Path(s3esti_path)
        _sys.path.insert(0, str(s3esti_path))
        from abso_esti_model.abso_esti_net import EstiNet
        from abso_esti_model.abso_esti_utils import get_wrap_patches_multi
        self._get_wrap_patches_multi = get_wrap_patches_multi

        ckpt_file = s3esti_path / "abso_esti_model" / checkpoint
        ckpt = torch.load(str(ckpt_file), map_location=self._device, weights_only=False)
        ckpt["model_scale"].pop("base", None)
        ckpt["model_angle"].pop("base", None)

        self._scale_net = EstiNet(
            need_bn=True, device=str(self._device), out_channels=300,
            patch_size=patch_size, scale_ratio=self._scale_ratios,
        )
        self._scale_net.load_state_dict(ckpt["model_scale"], strict=True)
        self._scale_net.eval().to(self._device)

        self._angle_net = EstiNet(
            need_bn=True, device=str(self._device), out_channels=360,
            patch_size=patch_size, scale_ratio=self._scale_ratios,
        )
        self._angle_net.load_state_dict(ckpt["model_angle"], strict=True)
        self._angle_net.eval().to(self._device)

    def _estimate_scale_angle(self, image_path, points):
        """Estimate scale and angle for keypoints in a single image."""
        import cv2
        import torch

        img_rgb = cv2.cvtColor(cv2.imread(str(image_path)), cv2.COLOR_BGR2RGB)
        img_t = (torch.from_numpy(img_rgb[np.newaxis].transpose(0, 3, 1, 2)).float() - 128) / 129
        img_t = img_t.to(self._device)

        pts = np.asarray(points, dtype=np.float32)
        sa_results, pid_list = [], []
        self._get_wrap_patches_multi(
            img_t, self._scale_net, self._angle_net, pts,
            self._num_patches, self._scale_ratios,
            scale_angle_result_list=sa_results, point_id_list=pid_list,
        )

        raw_scales = 1.0 / (sa_results[0] + 1e-16)
        raw_angles = -sa_results[1]
        pids = pid_list[0]

        # Take the first (best) estimate per point
        scales = np.zeros(len(pts))
        angles = np.zeros(len(pts))
        seen = set()
        for i, pid in enumerate(pids):
            if pid not in seen:
                seen.add(pid)
                scales[pid] = raw_scales[i]
                angles[pid] = raw_angles[i]
        return scales, angles

    def upgrade(self, result, *, path1, path2):
        scales1, angles1 = self._estimate_scale_angle(path1, result.points1)
        scales2, angles2 = self._estimate_scale_angle(path2, result.points2)

        return ScaleOrientationMatchResult(
            points1=result.points1, points2=result.points2, scores=result.scores,
            angles1=angles1, angles2=angles2,
            scales1=scales1, scales2=scales2,
        )


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

FEATURE_TYPES = {
    "sift": SIFT,
    "orb": ORB,
    "superpoint": SuperPoint,
    "disk": DISK,
    "aliked": ALIKED,
    "loma": LoMa,
    "affnet": AffNet,
    "denseaffine": DenseAffine,
    "pop": POP,
}


# ---------------------------------------------------------------------------
# CLI helpers
# ---------------------------------------------------------------------------

def add_feature_args(parser):
    """Add --features, --matcher, --num_features, --lowe_ratio, --resize to an ArgumentParser.
    Choices are derived from ``FeatureRegistry``, populated by ``__init_subclass__``."""
    parser.add_argument("--features", choices=FeatureRegistry.feature_names(), default="sift")
    parser.add_argument("--matcher", choices=FeatureRegistry.matcher_names(), default="bfmatcher",
                        help="OpenCV BFMatcher with Lowe's ratio test, or LightGlue")
    parser.add_argument("--num_features", type=int, default=None,
                        help="Keypoints per image (default: 4000 for BFMatcher, 2048 for LightGlue and LoMa)")
    parser.add_argument("--lowe_ratio", type=float, default=None,
                        help="Lowe's ratio test threshold for BFMatcher (default: 0.8)")
    parser.add_argument("--resize", type=int, default=None,
                        help="With LightGlue: resize the shorter side of the images to this size")


def create_matcher(args):
    """Create a FeatureMatcher from parsed CLI args.  Looks up the class from
    ``FeatureRegistry`` and calls its ``from_args`` classmethod.  Raises ValueError on
    invalid combinations."""
    if args.resize and args.matcher != "lightglue":
        raise ValueError("--resize needs --matcher lightglue")
    if args.lowe_ratio is not None and args.matcher != "bfmatcher":
        raise ValueError("--lowe_ratio needs --matcher bfmatcher")
    cls = FeatureRegistry.get(args.features)
    if cls is None:
        raise ValueError(f"unknown --features {args.features}")
    matchers = getattr(cls, 'cli_matchers', set())
    if matchers and args.matcher not in matchers:
        raise ValueError(f"--features {args.features} doesn't support --matcher {args.matcher}")
    return cls.from_args(args)
