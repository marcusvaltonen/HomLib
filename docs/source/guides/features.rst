.. _feature_types_guide:

Feature Types Guide
===================

HomLib supports three kinds of local features, each enabling different
solvers with different minimal sample sizes.  Richer features encode more
geometric information per correspondence, which reduces the minimal sample
size and dramatically improves RANSAC convergence — especially at high
outlier ratios.


Overview
--------

.. list-table::
   :header-rows: 1
   :widths: 18 25 25 12 20

   * - Feature type
     - Input data
     - Typical detectors
     - Min corrs.
     - Available solvers
   * - :attr:`~homlib.FeatureType.Points`
     - ``x``, ``y`` (2×N)
     - SIFT, ORB, SuperPoint
     - 5–6
     - 11 solvers
   * - :attr:`~homlib.FeatureType.Affine`
     - ``x``, ``y`` (2×N) + ``A`` (N×2×2)
     - Hessian-Affine, AffNet, MSER
     - 2
     - 2 solvers
   * - :attr:`~homlib.FeatureType.Orientation`
     - ``x``, ``y`` (2×N) + ``ori`` (2×N)
     - SIFT, ORB (with angle)
     - 3
     - 1 solver

A lower minimal sample size means RANSAC needs exponentially fewer iterations
for the same inlier ratio.  For example, at 50% inliers, a 2-point solver
needs roughly 4× the iterations of a single attempt, while a 5-point solver
needs roughly 32×.


Point features
--------------

The most widely supported feature type.  Any standard keypoint detector
and matcher produces point correspondences.

**Data format:**

- ``x`` — shape ``(2, N)``, centred normalised coordinates in image 1.
- ``y`` — shape ``(2, N)``, centred normalised coordinates in image 2.

**Centring and normalisation:**

.. code-block:: python

   import numpy as np

   # Raw pixel matches from a detector
   # pts1, pts2: shape (N, 2) in pixel coordinates
   cx, cy = width / 2, height / 2
   f = focal_length_pixels

   x = ((pts1 - [cx, cy]) / f).T   # shape (2, N)
   y = ((pts2 - [cx, cy]) / f).T

**Usage:**

.. code-block:: python

   import homlib

   solver = homlib.Wadenback3DV2026OneSided()
   poses = solver.estimate(x, y)

   # With RANSAC
   lomsac = homlib.LOMSAC(solver, inlier_threshold=0.07)
   estimate, stats = lomsac.estimate(x, y)


Affine-covariant features
--------------------------

Affine-covariant features capture the local affine deformation around each
keypoint.  Detectors like Hessian-Affine, MSER, or learned methods like
AffNet produce a 2×2 matrix describing the local image transformation at
each correspondence.

Because each match carries more geometric information, only **2
correspondences** are needed for a minimal sample.

**Data format:**

- ``x``, ``y`` — shape ``(2, N)``, same as point features.
- ``A`` — shape ``(N, 2, 2)``, local affine transformation matrices.
  ``A[i]`` is the 2×2 matrix relating the local patch in image 1 to image 2
  for correspondence ``i``.

**Constructing A from OpenCV keypoints:**

When using detectors that output elliptical regions (e.g. Hessian-Affine),
the affine matrix ``A`` encodes the transformation between the local coordinate
frames of matched features.

.. code-block:: python

   import numpy as np
   import homlib

   # A_matrices: shape (N, 2, 2) from your detector/matcher
   solver = homlib.ValtonenOrnhagICPR2026AffineOneSided()
   poses = solver.estimate(x, y, A=A_matrices)

   # With RANSAC
   lomsac = homlib.LOMSAC(solver, inlier_threshold=0.07)
   estimate, stats = lomsac.estimate(x, y, A=A_matrices)

**Available solvers:**

- :class:`~homlib.BarathVISAPP2016Affine` — 2-point, pure homography (no distortion).
- :class:`~homlib.ValtonenOrnhagICPR2026AffineOneSided` — 2-point, one-sided distortion (left).
- :class:`~homlib.ValtonenOrnhagICPR2026AffineOneSidedRight` — 2-point, one-sided distortion (right).


Orientation-covariant features
-------------------------------

Many keypoint detectors (SIFT, ORB, etc.) assign a dominant orientation
angle to each keypoint.  Orientation-covariant solvers exploit this
additional information, needing **3 correspondences** instead of 5.

**Data format:**

- ``x``, ``y`` — shape ``(2, N)``, same as point features.
- ``ori`` — shape ``(2, N)``, unit orientation vectors.  Column ``i`` is
  ``[cos(θ₁ᵢ - θ₂ᵢ), sin(θ₁ᵢ - θ₂ᵢ)]`` where ``θ₁ᵢ`` and ``θ₂ᵢ``
  are the keypoint angles in image 1 and image 2 respectively.

**Converting from SIFT/ORB angles:**

.. code-block:: python

   import numpy as np
   import homlib

   # angles1, angles2: shape (N,) in radians, from SIFT or ORB
   angle_diff = angles1 - angles2
   ori = np.vstack([np.cos(angle_diff), np.sin(angle_diff)])  # shape (2, N)

   solver = homlib.ValtonenOrnhagICPR2026OrientationOneSided()
   poses = solver.estimate(x, y, ori=ori)

   # With RANSAC
   lomsac = homlib.LOMSAC(solver, inlier_threshold=0.07)
   estimate, stats = lomsac.estimate(x, y, ori=ori)

**Available solvers:**

- :class:`~homlib.ValtonenOrnhagICPR2026OrientationOneSided` — 3-point,
  one-sided distortion (left).
- :class:`~homlib.ValtonenOrnhagICPR2026OrientationOneSidedRight` — 3-point,
  one-sided distortion (right).


Using features with RANSAC backends
------------------------------------

All three feature types work with all four RANSAC backends.  Pass the extra
feature data as keyword arguments to ``estimate()``:

.. code-block:: python

   import homlib

   # Point solver — no extra data needed
   solver_pt = homlib.Wadenback3DV2026OneSided()
   est, stats = homlib.LOMSAC(solver_pt, inlier_threshold=0.07).estimate(x, y)

   # Affine solver — pass A
   solver_aff = homlib.ValtonenOrnhagICPR2026AffineOneSided()
   est, stats = homlib.LOMSAC(solver_aff, inlier_threshold=0.07).estimate(x, y, A=A)

   # Orientation solver — pass ori
   solver_ori = homlib.ValtonenOrnhagICPR2026OrientationOneSided()
   est, stats = homlib.LOMSAC(solver_ori, inlier_threshold=0.07).estimate(x, y, ori=ori)

The same pattern works for :class:`~homlib.GCRANSAC`, :class:`~homlib.MAGSAC`,
and :class:`~homlib.SupeRANSAC`.

.. seealso::

   :ref:`solver_guide` for the full solver catalog, :ref:`backend_guide` for
   backend configuration, and the :ref:`examples` for worked examples.
