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
with the one-parameter **division model** :cite:`Fitzgibbon2001`:

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
   * - **One-sided (left)**
     - Distortion on the second image (**y**); the first image (**x**) is
       distortion-free. ``OneSided`` (no suffix) is the default.
     - ``distortion_parameter2``
   * - **One-sided (right)**
     - Distortion on the first image (**x**); the second image (**y**) is
       distortion-free. Solver names end in ``OneSidedRight``.
     - ``distortion_parameter``
   * - **Two-sided equal**
     - Same distortion coefficient on both images (:math:`k_1 = k_2`).
     - Both set to the same value
   * - **Two-sided independent**
     - Independent distortion on each image (:math:`k_1 \neq k_2`).
     - Both set independently
   * - **No distortion**
     - Pure homography, :math:`k = 0`.
     - Neither set

**Left vs. right** refers to which image carries the radial distortion in a
one-sided model.  "Left" (the default, ``OneSided``) places distortion on the
**y** argument (second image); "right" (``OneSidedRight``) places it on the
**x** argument (first image).  The ``OneSidedRight`` variants are derived
automatically from the corresponding ``OneSided`` solver, so they share the
same numerical core.  Choose whichever matches which of your two images is
distorted.


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
     - :cite:`Cai2025`
   * - :class:`~homlib.DLT`
     - 4+
     - None
     - :cite:`HartleyZisserman2003`
   * - :class:`~homlib.FitzgibbonCVPR2001OneSided`
     - 5+
     - One-sided (left)
     - :cite:`Fitzgibbon2001`
   * - :class:`~homlib.FitzgibbonCVPR2001OneSidedRight`
     - 5+
     - One-sided (right)
     - :cite:`Fitzgibbon2001`
   * - :class:`~homlib.FitzgibbonCVPR2001TwoSidedEqual`
     - 5+
     - Two-sided equal
     - :cite:`Fitzgibbon2001`
   * - :class:`~homlib.KukelovaCVPR2015TwoSidedEqual`
     - 5
     - Two-sided equal
     - :cite:`Kukelova2015`
   * - :class:`~homlib.KukelovaCVPR2015TwoSidedEqual6Pt`
     - 6
     - Two-sided equal
     - :cite:`Kukelova2015`
   * - :class:`~homlib.KukelovaCVPR2015TwoSided`
     - 5
     - Two-sided indep.
     - :cite:`Kukelova2015`
   * - :class:`~homlib.KukelovaCVPR2015TwoSided6Pt`
     - 6
     - Two-sided indep.
     - :cite:`Kukelova2015`
   * - :class:`~homlib.NakanoICPR2025OneSided`
     - 5+
     - One-sided (left)
     - :cite:`Nakano2024`
   * - :class:`~homlib.NakanoICPR2025OneSidedRight`
     - 5+
     - One-sided (right)
     - :cite:`Nakano2024`
   * - :class:`~homlib.Wadenback3DV2026OneSided`
     - 5
     - One-sided (left)
     - :cite:`Wadenback2026`
   * - :class:`~homlib.Wadenback3DV2026OneSidedRight`
     - 5
     - One-sided (right)
     - :cite:`Wadenback2026`
   * - :class:`~homlib.Wadenback3DV2026TwoSidedEqual`
     - 5
     - Two-sided equal
     - :cite:`Wadenback2026`
   * - :class:`~homlib.Wadenback3DV2026TwoSided`
     - 5
     - Two-sided indep.
     - :cite:`Wadenback2026`
   * - :class:`~homlib.GuoPRL2019`
     - 4
     - None
     - :cite:`Guo2019`

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
     - 2+
     - None
     - :cite:`Barath2016`
   * - :class:`~homlib.ValtonenOrnhagICPR2026AffineOneSided`
     - 2+
     - One-sided (left)
     - :cite:`ValtonenOrnhagAdalbjornsson2026`
   * - :class:`~homlib.ValtonenOrnhagICPR2026AffineOneSidedRight`
     - 2+
     - One-sided (right)
     - :cite:`ValtonenOrnhagAdalbjornsson2026`

Orientation-covariant solvers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 45 8 15 32

   * - Solver
     - Corrs.
     - Distortion
     - Reference
   * - :class:`~homlib.ValtonenOrnhagICPR2026OrientationOneSided`
     - 3+
     - One-sided (left)
     - :cite:`ValtonenOrnhagAdalbjornsson2026`
   * - :class:`~homlib.ValtonenOrnhagICPR2026OrientationOneSidedRight`
     - 3+
     - One-sided (right)
     - :cite:`ValtonenOrnhagAdalbjornsson2026`


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
   - Affine, one-sided: :class:`~homlib.ValtonenOrnhagICPR2026AffineOneSided`
   - Orientation, one-sided: :class:`~homlib.ValtonenOrnhagICPR2026OrientationOneSided`

   The Wadenbäck 3DV 2026 solvers are the newest and are generally the most
   numerically stable.  They enable ``extra_check`` by default, which filters
   spurious solutions.

Many of the solvers were generated using the automatic generator by
:cite:`Larsson2017`.

.. seealso::

   :ref:`api` for the full API reference, :ref:`references` for the original
   papers, and the :ref:`examples` for worked examples.


Sample check
------------

The 4-point solvers (:class:`~homlib.DLT`, :class:`~homlib.GuoPRL2019`,
:class:`~homlib.CaiPAMI2025`) can check their four points before computing a homography, and
return no solution for an invalid sample (RANSAC then moves on to the next one). The check is on
by default, in Python and in the robust estimators; in C++, ``get<false>`` skips it, which is
faster, but an invalid sample can then give a wrong homography. The Wadenbäck 3DV 2026 solvers
always check the undistorted points of each root, and drop the roots that give an invalid sample
(from five points, they first choose the four that are farthest from being collinear).

A sample is valid when:

- **no three points are collinear**, in either image: a triangle counts as collinear when its
  area is below :math:`10^{-10}` times the sum of the six squared distances between the four
  points. This only catches collinearity up to rounding errors (in particular four collinear
  points), and doesn't depend on the scale of the coordinates;
- **the four triangles all keep their orientation, or all flip it**, from one image to the other:
  otherwise the homography would map some of the points from behind a camera, which a plane seen
  by two cameras can't do.

Range check on the distortion
-----------------------------

Most solvers with distortion find it as the roots of a polynomial and return several
solutions, most with implausible distortion parameters. A range check keeps only the solutions
with both parameters in ``[k_min, k_max]`` (a parameter equal to 0, the one that a one-sided
solver doesn't estimate, passes). The solvers apply it to each root before computing
its homography, so RANSAC has fewer hypotheses to score. It only selects the hypotheses: the
refinement isn't constrained.

Where it applies
^^^^^^^^^^^^^^^^

The robust estimators (LOMSAC and the other backends) apply ``[-1, 0.1]`` by default. Before
calling a solver, they normalize each minimal sample: they compute
:math:`f_0 = \mathrm{mean}(|x|) / \sqrt{2}` from the first image's points :math:`x`, measured
from the distortion center (the origin), and divide the points of both images by :math:`f_0`
(the solvers with equal distortion in both images need the same scale in both). The range
applies in these normalized coordinates, whatever the units of the input.

:meth:`~homlib.Solver.set_k_range` replaces the range, ``set_k_range()`` without bounds
disables it, and :meth:`~homlib.Solver.reset_k_range` restores the default. The direct
:meth:`~homlib.Solver.estimate` calls the solver on the points as given, without a range unless
one is set; a range set with :meth:`~homlib.Solver.set_k_range` is then converted to the points'
scale (:math:`f_0` from the first image's points), or used as is with
``normalize_image_coord = False``.

Why the default is safe
^^^^^^^^^^^^^^^^^^^^^^^

Lens reviews give a lens's distortion as the relative displacement :math:`D` of the image
corners. To first order, :math:`D \approx k\, r_c^2`, where :math:`r_c` is the corners' radius in
the coordinates used. For points spread over an image, :math:`f_0` is about 0.19 times the image
diagonal (for any common aspect ratio), so in HomLib's normalized coordinates
:math:`r_c \approx 2.6` and :math:`k \approx D / 7`, whatever the focal length. ``[-1, 0.1]``
thus admits up to about 70% of pincushion distortion (:math:`k > 0`) at the corners, and any
barrel distortion (:math:`k < 0`) that the division model can represent there, since that requires
:math:`1 + k r_c^2 > 0`, i.e. :math:`k > -0.15`. Real lenses are far inside:

- **Smartphones and action cameras**: on more than 300 such cameras, López-Antequera et al.
  :cite:`LopezAntequera2019` (Fig. 3) measured the first coefficient of Brown's model in about
  ``[-0.4, 0]`` in *focal-normalized coordinates* (the centered coordinates divided by the focal
  length :math:`f`); to first order, the division model's parameter :math:`k_f` is the same. In the
  normalized coordinates, it is :math:`k = k_f (f_0 / f)^2`, where :math:`f_0`, about 0.19 times
  the image diagonal, equals the focal length of an 8 mm lens in 35 mm equivalent (43.3 mm
  diagonal): for a 14 mm action camera with :math:`k_f \approx -0.4`, :math:`k \approx -0.14`.
- **Long lenses**, which have pincushion distortion: lens reviews measured 0.8% at 200 mm for a
  70–200 mm lens and about 3% at the long end of a 20–200 mm superzoom [OpticalLimits], i.e.
  :math:`k \approx +0.001` to :math:`+0.004`.

The default thus only removes roots that no real lens would give. See the example
:ref:`sphx_glr_gallery_example_k_range.py`.

Turning off the normalization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With ``solver.normalize_image_coord = False`` (:attr:`~homlib.Solver.normalize_image_coord`),
the estimators use the points as given, and the range applies in their coordinates. The default
``[-1, 0.1]`` then applies in these units too: for pixel coordinates, it rejects almost nothing,
so set a range suited to them with :meth:`~homlib.Solver.set_k_range`.

This is useful when each image's focal length is known: divide its centered coordinates by it,
and the range is in focal-normalized coordinates. Choose it for the lens there: for the same
corner distortion, :math:`k_f` grows with the square of the focal length
(:math:`D \approx k_f (21.6\,\mathrm{mm} / f)^2` for a 35 mm-equivalent :math:`f`), so 1% of
pincushion distortion at the corners of a 200 mm lens is :math:`k_f \approx +0.9`, and 2% at
50 mm already :math:`k_f \approx +0.11`. ``[-1, 0.1]`` suits wide-angle lenses (barrel
distortion), not lenses with pincushion distortion. Dividing by the focal length has two uses:

- **A tighter range.** For wide-angle and action cameras (barrel distortion), a range such as
  ``[-0.6, 0.05]`` applies to all of them. In the setup of
  ``python/tests/test_normalization.py`` (a 14 mm lens, Wadenbäck et al.'s two-sided solver), it
  leaves about a third as many solutions per minimal sample as the default, and keeps the true
  solution.
- **Images with very different pixel scales.** Since the estimators take :math:`f_0` from the
  first image's points, a very large first image can put the second image's :math:`k` at the
  wrong scale. For a panorama with a focal length of 16000 px and a 1000 × 1000 view of it with a
  focal length of 1000 px, the view's true :math:`k` is outside ``[-1, 0.1]`` in every minimal sample, while
  with each image divided by its focal length it is inside.

HomLib uses one scale for both images because the solvers with equal distortion in both images
need it. If each image were divided by its own :math:`f_0`, the second case couldn't happen: a
barrel distortion would always satisfy :math:`|k| < 1/\max r^2 \le 1/2` there (since
:math:`1 + k r^2 > 0` at every distorted radius :math:`r`, and the mean radius is :math:`\sqrt{2}`).
