.. _solver_guide:

Solver Guide
============

HomLib provides 17 solvers from 9 papers for estimating homographies
with simultaneous radial distortion correction. This guide explains the
mathematical model, the available solvers, and how to choose the right one.


The radial distortion homography model
---------------------------------------

A **homography** :math:`\mathbf{H}` is a :math:`3 \times 3` matrix that maps
points between two images of a planar scene:

.. math::

   \tilde{\mathbf{y}} \sim \mathbf{H}\, \tilde{\mathbf{x}}

where :math:`\tilde{\mathbf{x}}` and :math:`\tilde{\mathbf{y}}` are homogeneous
coordinates.

**Radial distortion** bends straight lines in the image.  HomLib models this
with the one-parameter **division model** [Fitzgibbon2001]_:

.. math::

   \begin{pmatrix} \mathbf{x}_d \\ 1 + k\|\mathbf{x}_d\|^2 \end{pmatrix}
   \sim
   \begin{pmatrix} \mathbf{x}_u \\ 1 \end{pmatrix}

Equivalently, the undistorted coordinates are obtained by dividing the observed
(distorted) coordinates by :math:`1 + k\|\mathbf{x}_d\|^2`:

.. math::

   \mathbf{x}_u = \frac{\mathbf{x}_d}{1 + k\|\mathbf{x}_d\|^2}

Here :math:`\mathbf{x}_d` are the *distorted* (observed) coordinates,
:math:`\mathbf{x}_u` are the *undistorted* (ideal) coordinates, and
:math:`k` is the distortion parameter.  Negative :math:`k` produces
barrel distortion; positive :math:`k` produces pincushion distortion.

HomLib solvers **simultaneously** estimate :math:`\mathbf{H}` and :math:`k`
from point correspondences, without requiring prior lens calibration.


Coordinate conventions
^^^^^^^^^^^^^^^^^^^^^^

Solvers that estimate radial distortion require **centered** coordinates::

   x_centered = (x_pixel - cx) / s
   y_centered = (y_pixel - cy) / s

where ``(cx, cy)`` is the principal point (typically the image centre) and
``s`` is an arbitrary normalization scalar (e.g. the image width, the focal
length, or any convenient value).  Camera calibration is **not** required —
any scalar that brings coordinates to a reasonable range works.  Centering
matters because the division model assumes the distortion centre is at the
origin; inlier thresholds are expressed in these normalised units.

Pure homography solvers (no distortion) do not require centering, but
consistent normalization still improves numerical conditioning.

The output :class:`~homlib.PoseData` contains:

- ``homography`` — the :math:`3 \times 3` homography matrix.
- ``distortion_parameter`` — :math:`k_1` (first image).
- ``distortion_parameter2`` — :math:`k_2` (second image).
- ``focal_length`` — currently unused; reserved for future solvers.


Distortion configurations
^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 50 25

   * - Configuration
     - Description
     - PoseData fields set
   * - **One-sided**
     - Distortion on one image only (the other is assumed distortion-free).
     - ``distortion_parameter`` or ``distortion_parameter2``
   * - **Two-sided equal**
     - Same distortion coefficient on both images (:math:`k_1 = k_2`).
     - Both set to the same value
   * - **Two-sided independent**
     - Independent distortion on each image (:math:`k_1 \neq k_2`).
     - Both set independently
   * - **No distortion**
     - Pure homography, :math:`k = 0`.
     - Neither set


.. _feature_types_overview:

Feature types
^^^^^^^^^^^^^

HomLib solvers accept three kinds of features, distinguished by the
:class:`~homlib.FeatureType` enum:

- :attr:`~homlib.FeatureType.Points` — standard keypoint matches ``(x, y)``.
  Requires 5–6 correspondences per minimal sample.
- :attr:`~homlib.FeatureType.Affine` — affine-covariant features including
  a local 2×2 transformation ``A``.  Requires only **2** correspondences.
- :attr:`~homlib.FeatureType.Orientation` — orientation-covariant features
  with a keypoint angle.  Requires **3** correspondences.

Richer features need fewer correspondences, which helps RANSAC converge
faster (exponentially fewer iterations for the same inlier ratio).
See :ref:`feature_types_guide` for a detailed explanation.


Solver catalog
--------------

Point-based solvers
^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 35 8 15 42

   * - Solver
     - Corrs.
     - Distortion
     - Reference
   * - :class:`~homlib.CaiPAMI2025`
     - 4
     - None
     - [Cai2025]_
   * - :class:`~homlib.DLT`
     - 4
     - None
     - [HartleyZisserman2003]_
   * - :class:`~homlib.FitzgibbonCVPR2001OneSided`
     - 5
     - One-sided
     - [Fitzgibbon2001]_
   * - :class:`~homlib.FitzgibbonCVPR2001TwoSidedEqual`
     - 5
     - Two-sided equal
     - [Fitzgibbon2001]_
   * - :class:`~homlib.KukelovaCVPR2015TwoSidedEqual`
     - 5
     - Two-sided equal
     - [Kukelova2015]_
   * - :class:`~homlib.KukelovaCVPR2015TwoSidedEqual6Pt`
     - 6
     - Two-sided equal
     - [Kukelova2015]_
   * - :class:`~homlib.KukelovaCVPR2015TwoSided`
     - 5
     - Two-sided indep.
     - [Kukelova2015]_
   * - :class:`~homlib.KukelovaCVPR2015TwoSided6Pt`
     - 6
     - Two-sided indep.
     - [Kukelova2015]_
   * - :class:`~homlib.NakanoICPR2025OneSided`
     - 5
     - One-sided (left)
     - [Nakano2024]_
   * - :class:`~homlib.NakanoICPR2025OneSidedRight`
     - 5
     - One-sided (right)
     - [Nakano2024]_
   * - :class:`~homlib.Wadenback3DV2026OneSided`
     - 5
     - One-sided
     - [Wadenback2026]_
   * - :class:`~homlib.Wadenback3DV2026TwoSidedEqual`
     - 5
     - Two-sided equal
     - [Wadenback2026]_
   * - :class:`~homlib.Wadenback3DV2026TwoSided`
     - 5
     - Two-sided indep.
     - [Wadenback2026]_
   * - :class:`~homlib.GuoPRL2019`
     - 4
     - None
     - [Guo2019]_

Affine-covariant solvers
^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 40 8 15 37

   * - Solver
     - Corrs.
     - Distortion
     - Reference
   * - :class:`~homlib.BarathVISAPP2016Affine`
     - 2
     - None
     - [Barath2016]_
   * - :class:`~homlib.ValtonenOrnhagICPR2026AffineSingleSided`
     - 2
     - One-sided
     - [ValtonenOrnhagAdalbjornsson2026]_

Orientation-covariant solvers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 45 8 15 32

   * - Solver
     - Corrs.
     - Distortion
     - Reference
   * - :class:`~homlib.ValtonenOrnhagICPR2026OrientationSingleSided`
     - 3
     - One-sided
     - [ValtonenOrnhagAdalbjornsson2026]_


Choosing a solver
-----------------

1. **What feature data do you have?**

   - Standard keypoint matches only → use a **point-based** solver.
   - Local affine frames (e.g. from AffNet, Hessian-Affine) → use an **affine** solver.
   - Keypoint orientations (e.g. from SIFT) → use an **orientation** solver.

2. **What distortion model do you need?**

   - One camera is distortion-free (e.g. rectified or pinhole) → **one-sided**.
   - Both cameras have the same lens → **two-sided equal**.
   - Different lenses on each camera → **two-sided independent**.
   - No distortion at all → :class:`~homlib.CaiPAMI2025` (point features) or
     :class:`~homlib.BarathVISAPP2016Affine` (affine features).

3. **Recommended defaults:**

   - Point, one-sided: :class:`~homlib.Wadenback3DV2026OneSided`
   - Point, two-sided equal: :class:`~homlib.Wadenback3DV2026TwoSidedEqual`
   - Point, two-sided independent: :class:`~homlib.Wadenback3DV2026TwoSided`
   - Affine, one-sided: :class:`~homlib.ValtonenOrnhagICPR2026AffineSingleSided`
   - Orientation, one-sided: :class:`~homlib.ValtonenOrnhagICPR2026OrientationSingleSided`

   The Wadenbäck 3DV 2026 solvers are the newest and are generally the most
   numerically stable.  They enable ``extra_check`` by default, which filters
   spurious solutions.

Many of the solvers were generated using the automatic generator by
[Larsson2017]_.

.. seealso::

   :ref:`api` for the full API reference, :ref:`references` for the original
   papers, and the :ref:`examples` for worked examples.
