.. _api:

homlib API
=================

.. currentmodule:: homlib

Architecture
------------

HomLib uses a **composable** design: solvers and RANSAC backends are independent.
Any solver can be paired with any RANSAC backend.

- **Solvers** estimate homographies (and radial distortion parameters) from a
  minimal point sample. Each solver is a class named after its paper, e.g.
  :class:`Wadenback3DV2026OneSided`. Calling :meth:`~Solver.estimate` returns a
  list of algebraic solutions.

- **RANSAC backends** wrap a solver and run a robust estimation loop in C++.
  :class:`LOMSAC` is always available; :class:`SupeRANSAC` is optional.
  Both expose a single :meth:`~RobustEstimator.estimate` method.

.. code-block:: python

   solver = homlib.Wadenback3DV2026OneSided()

   # Direct (non-robust) estimation
   poses = solver.estimate(x, y)

   # Robust with LOMSAC
   lomsac = homlib.LOMSAC(solver, squared_inlier_threshold=0.005)
   estimate, stats = lomsac.estimate(x, y)

   # Robust with SupeRANSAC
   sa = homlib.SupeRANSAC(solver, inlier_threshold=0.07)
   estimate, stats = sa.estimate(x, y)

Data types
----------

.. autosummary::
   :toctree: generated/

    PoseData
    LORansacOptions
    RansacStatistics
    DistortionCase
    ProblemConfig
    ProblemInstance

Enums
-----

.. autosummary::
   :toctree: generated/

    SolverId
    FeatureType

Base classes
------------

.. autoclass:: Solver
   :members:
   :undoc-members:

.. autoclass:: RobustEstimator
   :members:

Point-based solvers
-------------------

Direct solver classes that return all algebraic solutions from a minimal sample.
Distortion suffixes: ``OneSided`` — single distortion coefficient on one image,
``TwoSidedEqual`` — equal distortion on both images, ``TwoSided`` — independent
distortion on each image.

All point-based solvers have ``feature_type == FeatureType.Points``.

.. autoclass:: FitzgibbonCVPR2001OneSided
.. autoclass:: FitzgibbonCVPR2001TwoSidedEqual
.. autoclass:: KukelovaCVPR2015TwoSidedEqual
.. autoclass:: KukelovaCVPR2015TwoSidedEqual6Pt
.. autoclass:: KukelovaCVPR2015TwoSided
.. autoclass:: KukelovaCVPR2015TwoSided6Pt
.. autoclass:: NakanoICPR2025OneSided
.. autoclass:: NakanoICPR2025OneSidedRight
.. autoclass:: Wadenback3DV2026OneSided
.. autoclass:: Wadenback3DV2026TwoSidedEqual
.. autoclass:: Wadenback3DV2026TwoSided

Affine-covariant solvers
------------------------

Solvers that leverage affine-covariant features, requiring fewer correspondences
than point-based solvers. These have ``feature_type == FeatureType.Affine`` and
require the ``A`` parameter in :meth:`~Solver.estimate`.

.. autoclass:: BarathVISAPP2016Affine
.. autoclass:: ValtonenOrnhagICPR2026AffineSingleSided

Orientation-covariant solvers
-----------------------------

Solvers that leverage orientation-covariant features. These have
``feature_type == FeatureType.Orientation`` and require the ``ori`` parameter
in :meth:`~Solver.estimate`.

.. autoclass:: ValtonenOrnhagICPR2026OrientationSingleSided

RANSAC backends
---------------

Robust estimation backends wrap a solver and run a RANSAC-style loop entirely
in C++. All backends inherit from :class:`RobustEstimator`.

LOMSAC
^^^^^^

.. autoclass:: LOMSAC
   :members:
   :undoc-members:

SupeRANSAC
^^^^^^^^^^

.. autoclass:: SupeRANSAC
   :members:
   :undoc-members:

.. note::

   :class:`SupeRANSAC` requires HomLib to be built with
   ``-DHOMLIB_USE_SUPERANSAC=ON``. See :ref:`dependencies` for required
   system libraries (Boost, OpenCV). If the backend is not available,
   constructing a :class:`SupeRANSAC` instance raises :exc:`ImportError`.

   Check availability at runtime with ``homlib._HAS_SUPERANSAC``.

SupeRANSAC configuration types
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These types are only available when HomLib is built with SupeRANSAC support.

- **SuperansacSettings** — full configuration object (inlier threshold, scoring
  type, sampler type, local optimisation, etc.). Pass as ``settings=`` to
  :class:`SupeRANSAC`, or set individual fields via keyword arguments.
- **SuperansacStatistics** — results returned alongside the best model
  (``num_iterations``, ``inliers``).
- **SuperansacScoringType** — ``MSAC``, ``MAGSAC``, etc.
- **SuperansacSamplerType** — ``Uniform``, ``PROSAC``, etc.
- **SuperansacLocalOptType** — ``None_``, ``IRLS``, ``NestedRANSAC``,
  ``LSQ``, ``GCRANSAC``, etc.

Synthetic data generation
-------------------------

Utilities for generating synthetic problem instances with known ground truth,
useful for benchmarking and testing.

.. autosummary::
   :toctree: generated/

    generate_problem_instance
