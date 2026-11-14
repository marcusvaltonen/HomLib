.. _feature_matchers:

Feature Matchers and Upgraders
==============================

The benchmark framework provides a unified interface for feature extraction,
matching, and enrichment through the :class:`FeatureMatcher` and
:class:`FeatureUpgrader` base classes (in ``python/benchmark/common.py``).
Concrete implementations live in ``python/benchmark/features.py``.


Feature matchers
----------------

Each matcher is named by the **feature type** it detects, not by the matching
backend.  All matchers implement :meth:`extract` and :meth:`_match_features`,
and the inherited :meth:`match(path1, path2)` orchestrates both steps.
End-to-end systems (LoMa) override :meth:`match` directly.

.. list-table::
   :header-rows: 1
   :widths: 14 30 20 18 18

   * - Matcher
     - Supported feature types
     - Matcher backend
     - Returns
     - Reference
   * - :ref:`SIFT <fm_sift>`
     - Points, Orientation, ScaleAndOrientation
     - BFMatcher + ratio test
     - ``ScaleOrientationMatchResult``
     - :cite:`Lowe2004`
   * - :ref:`ORB <fm_orb>`
     - Points, Orientation
     - BFMatcher (Hamming) + ratio test
     - ``OrientationMatchResult``
     - :cite:`Rublee2011`
   * - :ref:`SuperPoint <fm_superpoint>`
     - Points
     - LightGlue / SuperGlue
     - ``PointMatchResult``
     - :cite:`DeTone2018`
   * - :ref:`DISK <fm_disk>`
     - Points
     - LightGlue
     - ``PointMatchResult``
     - :cite:`Tyszkiewicz2020`
   * - :ref:`ALIKED <fm_aliked>`
     - Points
     - LightGlue
     - ``PointMatchResult``
     - :cite:`Zhao2023`
   * - :ref:`LoMa <fm_loma>`
     - Points
     - End-to-end
     - ``PointMatchResult``
     - :cite:`Nordstrom2026`
   * - :ref:`AffNet <fm_affnet>`
     - Points, Orientation, Affine
     - BFMatcher + ratio test
     - ``AffineMatchResult``
     - :cite:`Mishkin2018`
   * - :ref:`DenseAffine <fm_denseaffine>`
     - Points, Orientation, Affine
     - End-to-end (stub)
     - ``AffineMatchResult``
     - :cite:`Stilcrad2024`
   * - :ref:`POP <fm_pop>`
     - Points
     - BFMatcher + ratio test (stub)
     - ``PointMatchResult``
     -


.. _fm_sift:

SIFT
^^^^

Scale-Invariant Feature Transform :cite:`Lowe2004`.

- **Supported feature types:** Points, Orientation, ScaleAndOrientation
- **Detector/descriptor:** OpenCV ``SIFT_create``
- **Matcher:** BFMatcher + Lowe's ratio test (default threshold 0.8)
- **Returns:** ``ScaleOrientationMatchResult`` with angles (from ``cv2.KeyPoint.angle``,
  converted to radians) and scales (from ``cv2.KeyPoint.size``)
- **Constructor:** ``SIFT(num_features=4000, ratio=0.8)``

With OpenCV 4.8+, uses SIFT's precise upscaling for improved keypoint accuracy.

.. code-block:: python

   from features import SIFT
   matcher = SIFT(num_features=4000, ratio=0.8)
   result = matcher("image1.jpg", "image2.jpg")
   # result.angles1, result.angles2: N arrays of angles in radians
   # result.scales1, result.scales2: N arrays of keypoint sizes


.. _fm_orb:

ORB
^^^

Oriented FAST and Rotated BRIEF :cite:`Rublee2011`.

- **Supported feature types:** Points, Orientation
- **Detector/descriptor:** OpenCV ``ORB_create``
- **Matcher:** BFMatcher with Hamming distance + ratio test
- **Returns:** ``OrientationMatchResult`` with angles (no meaningful scale)
- **Constructor:** ``ORB(num_features=4000, ratio=0.8)``


.. _fm_superpoint:

SuperPoint
^^^^^^^^^^

Self-supervised interest point detector and descriptor :cite:`DeTone2018`,
matched with LightGlue :cite:`Lindenberger2023` (default) or SuperGlue.

- **Supported feature types:** Points
- **Matcher:** LightGlue (default) or SuperGlue, selectable via ``matcher`` parameter
- **Returns:** ``PointMatchResult`` with match confidence scores
- **Constructor:** ``SuperPoint(num_features=2048, resize=None, matcher="lightglue")``

.. note::
   SuperPoint's weights are for non-commercial use only.
   Requires `LightGlue <https://github.com/cvg/LightGlue>`_.

.. code-block:: python

   from features import SuperPoint
   matcher = SuperPoint(num_features=1024, resize=480)
   result = matcher("image1.jpg", "image2.jpg")


.. _fm_disk:

DISK
^^^^

Learned local features with policy gradient :cite:`Tyszkiewicz2020`,
matched with LightGlue :cite:`Lindenberger2023`.

- **Supported feature types:** Points
- **Returns:** ``PointMatchResult``
- **Constructor:** ``DISK(num_features=2048, resize=None)``


.. _fm_aliked:

ALIKED
^^^^^^

Lighter keypoint and descriptor extraction via deformable transformation
:cite:`Zhao2023`, matched with LightGlue :cite:`Lindenberger2023`.

- **Supported feature types:** Points
- **Returns:** ``PointMatchResult``
- **Constructor:** ``ALIKED(num_features=2048, resize=None)``


.. _fm_loma:

LoMa
^^^^

Lossless Matching :cite:`Nordstrom2026` — an end-to-end system that detects,
describes, and matches in one forward pass.

- **Supported feature types:** Points
- **Returns:** ``PointMatchResult``
- **Constructor:** ``LoMa(num_keypoints=2048)``
- **Note:** Overrides ``match()`` directly (no separate extract/match steps).
  Pixel origin is the top-left corner: the framework shifts points by -0.5 px
  to OpenCV's convention.

Requires `LoMa <https://github.com/davnords/LoMa>`_ (``pip install -e .``
in a clone; the PyPI package named ``loma`` is unrelated).


.. _fm_affnet:

AffNet + HardNet
^^^^^^^^^^^^^^^^

AffNet :cite:`Mishkin2018` estimates an affine shape for each SIFT keypoint,
and HardNet :cite:`Mishchuk2017` computes descriptors on the affine-normalised
patches.  Matched with BFMatcher + ratio test.

- **Supported feature types:** Points, Orientation, Affine
- **Detector:** ``OpenCVDetectorWithAffNetKornia`` (kornia) on top of SIFT
- **Descriptor:** ``LAFDescriptor(HardNet)`` (kornia)
- **Returns:** ``AffineMatchResult`` (points, angles, scales, and full 2x2 affine frames)
- **Constructor:** ``AffNet(num_features=4000, ratio=0.8)``

Requires ``kornia`` and ``kornia_moons``.

.. code-block:: python

   from features import AffNet
   matcher = AffNet(num_features=4000)
   result = matcher("image1.jpg", "image2.jpg")
   A_relative = result.relative_affines()  # Nx2x2 for estimator.estimate(x, y, A=...)


.. _fm_denseaffine:

DenseAffine
^^^^^^^^^^^

Dense affine matching pipeline :cite:`Stilcrad2024` (DKMv3 + EstiNet + HardNet).

- **Supported feature types:** Points, Orientation, Affine
- **Status:** Stub — requires the DenseAffine repository.
- **Constructor:** ``DenseAffine(repo_path, weights_path=None, device=None)``


.. _fm_pop:

POP
^^^

Learned unsupervised interest point detector/descriptor.

- **Supported feature types:** Points
- **Status:** Stub — requires the POP repository.
- **Constructor:** ``POP(model_path=None, num_features=1000, device=None)``


Match result hierarchy
----------------------

All matchers return a subclass of ``PointMatchResult``.  The hierarchy encodes
which geometric fields are available:

.. code-block:: text

   PointMatchResult
    └── OrientationMatchResult         (+angles1, angles2)
         └── ScaleOrientationMatchResult   (+scales1, scales2)
              └── AffineMatchResult            (+affines1, affines2)

Each level inherits from the one above, so an ``AffineMatchResult`` can be used
anywhere a ``PointMatchResult`` is expected.  The ``AffineMatchResult`` class
also provides:

- ``relative_affines()`` — computes ``affines2 @ inv(affines1)``, the N×2×2
  input for ``estimator.estimate(x, y, A=...)``.
- ``from_lafs(laf1, laf2, scores=None)`` — class method to extract all fields
  from N×2×3 Local Affine Frames (LAFs).


Feature upgraders
-----------------

Upgraders enrich existing match results with additional geometric information.
They implement the ``FeatureUpgrader`` interface with a single ``upgrade()``
method.


AffineDecompositionUpgrader
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Downgrades ``AffineMatchResult`` to ``ScaleOrientationMatchResult`` by
extracting angles and scales from the individual affine frames (dropping the
full affine data).  Useful when you want to use an orientation or
scale+orientation solver with data that was originally matched as affine.

- **Input:** ``AffineMatchResult``
- **Output:** ``ScaleOrientationMatchResult``
- **No extra arguments needed.**


BarathAffineUpgrader
^^^^^^^^^^^^^^^^^^^^

Upgrades ``ScaleOrientationMatchResult`` to ``AffineMatchResult`` using the
epipolar constraint :cite:`Barath2023affine`.  Given known SIFT angles and
scales plus a fundamental matrix *F*, solves for the unknown affine parameters
(scale along *u*, scale along *v*, shear) in closed form per correspondence
(Eq. 11--18 of the paper).

The affine is decomposed as :math:`A = R_2 \cdot U \cdot R_1^T` where
:math:`U = \begin{bmatrix} q_u & w \\ 0 & q_v \end{bmatrix}` and
:math:`q_u = q_2^2 / (q_1^2 \cdot q_v)`.

- **Input:** ``ScaleOrientationMatchResult``
- **Output:** ``AffineMatchResult``
- **Required kwarg:** ``F`` — the 3×3 fundamental matrix.

.. code-block:: python

   from features import BarathAffineUpgrader
   upgrader = BarathAffineUpgrader()
   affine_result = upgrader.upgrade(sift_result, F=F_matrix)


S3EstiUpgrader
^^^^^^^^^^^^^^

Upgrades ``PointMatchResult`` to ``ScaleOrientationMatchResult`` using EstiNet,
a CNN that estimates per-keypoint scale (300 bins) and orientation (360 bins)
from multi-scale image patches.

- **Input:** ``PointMatchResult``
- **Output:** ``ScaleOrientationMatchResult``
- **Required kwargs:** ``path1``, ``path2`` — image paths for patch extraction.
- **Constructor:** ``S3EstiUpgrader(s3esti_path, checkpoint="S3Esti_ep30.pth", device=None, patch_size=32, scale_ratios=(0.5, 1, 2), num_patches=4)``

Requires the S3Esti repository and model weights.

.. code-block:: python

   from features import S3EstiUpgrader
   upgrader = S3EstiUpgrader("thirdparty/S3Esti")
   result = upgrader.upgrade(point_result, path1="img1.jpg", path2="img2.jpg")


Feature pipeline
----------------

Matchers and upgraders compose into pipelines:

.. code-block:: python

   from features import SIFT, BarathAffineUpgrader

   # Extract SIFT features (scale + orientation)
   sift = SIFT(num_features=4000)
   result = sift("image1.jpg", "image2.jpg")

   # Upgrade to affine using fundamental matrix
   upgrader = BarathAffineUpgrader()
   affine_result = upgrader.upgrade(result, F=F_matrix)

   # Use with an affine solver
   import homlib
   solver = homlib.ValtonenOrnhagICPR2026AffineOneSided()
   lomsac = homlib.LOMSAC(solver, inlier_threshold=0.07)
   x = affine_result.points1.T  # (2, N)
   y = affine_result.points2.T
   A = affine_result.relative_affines()
   estimate, stats = lomsac.estimate(x, y, A=A)
