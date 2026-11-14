"""Matcher classes: given two sets of extracted features, find correspondences.

Each matcher returns (idx1, idx2, scores): arrays of matched indices into the
two feature sets and a per-match score (higher = better).
"""

from abc import ABC, abstractmethod

import numpy as np


class Matcher(ABC):
    """A matching strategy that pairs features from two images."""

    @abstractmethod
    def match(self, features1, features2):
        """Returns (idx1, idx2, scores) as numpy arrays."""
        ...


class BFMatcher(Matcher):
    """OpenCV BFMatcher with Lowe's ratio test."""

    def __init__(self, lowe_ratio=0.8, norm=None):
        import cv2
        self._lowe_ratio = lowe_ratio
        self._bf = cv2.BFMatcher() if norm is None else cv2.BFMatcher(norm)

    def match(self, features1, features2):
        d1, d2 = features1["descriptors"], features2["descriptors"]
        if d1 is None or d2 is None or len(d1) < 2 or len(d2) < 2:
            return np.array([], dtype=int), np.array([], dtype=int), np.array([])
        matches = self._bf.knnMatch(d1, d2, k=2)
        good = [(m, n) for m, n in matches if m.distance < self._lowe_ratio * n.distance]
        if not good:
            return np.array([], dtype=int), np.array([], dtype=int), np.array([])
        idx1 = np.array([m.queryIdx for m, _ in good])
        idx2 = np.array([m.trainIdx for m, _ in good])
        scores = np.array([1 - m.distance / n.distance for m, n in good])
        return idx1, idx2, scores


class LightGlueMatcher(Matcher):
    """LightGlue matcher (Lindenberger et al., ICCV 2023)."""

    def __init__(self, feature_name, device):
        from lightglue import LightGlue
        self._matcher = LightGlue(features=feature_name).eval().to(device)

    def match(self, features1, features2):
        import torch
        from lightglue.utils import rbd
        with torch.no_grad():
            out = rbd(self._matcher({"image0": features1, "image1": features2}))
        m = out["matches"]
        return (m[:, 0].cpu().numpy(), m[:, 1].cpu().numpy(),
                out["scores"].cpu().numpy().astype(float))
