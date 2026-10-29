Welcome to homlib's documentation!
=========================================================

**HomLib** is a C++ library (with Python bindings) for computing homographies with
simultaneous radial distortion correction, optionally incorporating IMU data.
It supports point-based, affine-covariant, and orientation-covariant features.

.. code-block:: python

   import numpy as np
   import homlib

   # Choose a solver (each follows the Author-Venue-Year naming convention)
   solver = homlib.Wadenback3DV2026OneSided()

   # Direct estimation from a minimal sample
   poses = solver.estimate(x, y)
   print(poses[0].homography)

   # Robust estimation with LOMSAC (always available)
   lomsac = homlib.LOMSAC(solver, squared_inlier_threshold=0.005)
   estimate, stats = lomsac.estimate(x_noisy, y_noisy)

   # Robust estimation with SupeRANSAC (optional backend)
   sa = homlib.SupeRANSAC(solver, inlier_threshold=0.07)
   estimate, stats = sa.estimate(x_noisy, y_noisy)

Install from PyPI:

.. code-block:: bash

   pip install homlib

Getting started
---------------

HomLib has a composable architecture built around **solvers** and **RANSAC backends**.

**Solvers** implement specific minimal-solver algorithms from the literature.
Each solver is a class named after the paper (e.g. ``Wadenback3DV2026OneSided``).
Solvers have a ``feature_type`` property (:class:`~homlib.FeatureType`) indicating
whether they use point, affine, or orientation features:

.. code-block:: python

   solver = homlib.Wadenback3DV2026OneSided()      # point-based
   aff = homlib.ValtonenOrnhagICPR2026AffineSingleSided()  # affine
   ori = homlib.ValtonenOrnhagICPR2026OrientationSingleSided()  # orientation

**RANSAC backends** wrap a solver and run a robust estimation loop entirely in
C++. Any solver can be paired with any backend:

- :class:`~homlib.LOMSAC` — Locally Optimized MSAC (always available).
- :class:`~homlib.SupeRANSAC` — state-of-the-art RANSAC with configurable
  sampler, scoring, and local optimisation. Requires building with
  ``-DHOMLIB_USE_SUPERANSAC=ON``.

.. code-block:: python

   # Same solver, different backends
   solver = homlib.FitzgibbonCVPR2001OneSided()

   lomsac = homlib.LOMSAC(solver, squared_inlier_threshold=0.005)
   est_l, stats_l = lomsac.estimate(x, y)

   sa = homlib.SupeRANSAC(solver, inlier_threshold=0.07)
   est_s, stats_s = sa.estimate(x, y)

For affine or orientation solvers, pass the extra features to ``estimate()``:

.. code-block:: python

   aff_solver = homlib.ValtonenOrnhagICPR2026AffineSingleSided()
   lomsac = homlib.LOMSAC(aff_solver, squared_inlier_threshold=0.005)
   estimate, stats = lomsac.estimate(x, y, A=affine_matrices)

| `GitHub <https://github.com/marcusvaltonen/HomLib>`__ · `PyPI <https://pypi.org/project/homlib/>`__ · :ref:`api`

**Maintainer**: `Marcus Valtonen Örnhag <marcusvaltonen@gmail.com>`_

Citing HomLib
-------------

If you use HomLib in your research, please cite the relevant paper(s) from
the :ref:`references` page. For the library itself:

.. code-block:: bibtex

   @misc{HomLib,
       title  = {HomLib -- Homography Library of state-of-the-art methods},
       author = {Marcus Valtonen {\"O}rnhag},
       URL    = {https://github.com/marcusvaltonen/HomLib},
       year   = {2020},
   }


.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Getting started:

   installation/dependencies.rst
   installation/installation_users.rst
   installation/installation_developers.rst

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Examples and usage:

   gallery/index.rst

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Reference documentation:

   api/index.rst
   references.rst
