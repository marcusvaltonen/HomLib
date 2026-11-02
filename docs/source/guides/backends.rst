.. _backend_guide:

Backend Guide
=============

HomLib's RANSAC backends wrap a solver and run a complete robust estimation loop
in C++. The architecture is **composable**: any solver can be paired with any
backend. All backends inherit from :class:`~homlib.RobustEstimator` and expose a
single ``estimate()`` method.

.. code-block:: python

   import homlib

   solver = homlib.Wadenback3DV2026OneSided()

   # Same solver, different backends
   lomsac = homlib.LOMSAC(solver, inlier_threshold=0.07)
   gcransac = homlib.GCRANSAC(solver, threshold=0.07)
   magsac = homlib.MAGSAC(solver, sigma_max=10.0)
   superansac = homlib.SupeRANSAC(solver, inlier_threshold=0.07)


Backend comparison
------------------

.. list-table::
   :header-rows: 1
   :widths: 14 12 16 20 20 18

   * - Backend
     - Available
     - Threshold
     - Scoring
     - Local optimisation
     - Reference
   * - :class:`~homlib.LOMSAC`
     - Always
     - Squared
     - MSAC
     - LO-RANSAC + GN refinement
     - [Lebeda2012]_
   * - :class:`~homlib.GCRANSAC`
     - Optional
     - Raw pixel
     - MSAC
     - Graph-cut + IRLS
     - [Barath2018]_
   * - :class:`~homlib.MAGSAC`
     - Optional
     - σ\ :sub:`max`
     - Sigma consensus
     - IRLS (weighted)
     - [Barath2020]_
   * - :class:`~homlib.SupeRANSAC`
     - Optional
     - Raw pixel
     - Configurable
     - Configurable
     - [Barath2025]_


LOMSAC
------

**Locally Optimised MSAC** — always available, no extra build dependencies.

LOMSAC uses a squared inlier threshold and MSAC scoring. After finding a good
model, it refines via iterative least-squares on the inlier set (LO-RANSAC).

.. code-block:: python

   solver = homlib.Wadenback3DV2026OneSided()
   lomsac = homlib.LOMSAC(solver, inlier_threshold=0.07)
   estimate, stats = lomsac.estimate(x, y)

**Key settings** (via :class:`~homlib.LORansacOptions`):

.. list-table::
   :header-rows: 1
   :widths: 35 15 50

   * - Parameter
     - Default
     - Description
   * - ``inlier_threshold``
     - —
     - Inlier threshold (convenience; internally squared)
   * - ``squared_inlier_threshold``
     - —
     - Squared distance threshold (alternative to ``inlier_threshold``)
   * - ``max_num_iterations``
     - 10 000
     - Maximum RANSAC iterations
   * - ``success_probability``
     - 0.9999
     - Early termination probability
   * - ``num_lo_steps``
     - 4
     - Local optimisation refinement steps
   * - ``final_least_squares``
     - ``True``
     - Final Gauss-Newton refinement on all inliers
   * - ``random_seed``
     - 0
     - Random seed (0 = non-deterministic)

**When to use:** Default choice for production. Deterministic with
``random_seed``. Well-tested and fast.


GC-RANSAC
---------

**Graph-Cut RANSAC** [Barath2018]_ — uses graph-cut spatial coherence
to ensure inliers form spatially consistent regions. This helps reject
outliers that are close in residual but spatially inconsistent.

.. code-block:: python

   solver = homlib.Wadenback3DV2026OneSided()
   gc = homlib.GCRANSAC(solver, threshold=0.07)
   estimate, stats = gc.estimate(x, y)

Requires building with ``-DHOMLIB_USE_GCRANSAC=ON``. Check availability:

.. code-block:: python

   if homlib._HAS_GCRANSAC:
       gc = homlib.GCRANSAC(solver, threshold=0.07)

**Key settings** (via :class:`~homlib.GCRANSACSettings`):

.. list-table::
   :header-rows: 1
   :widths: 35 15 50

   * - Parameter
     - Default
     - Description
   * - ``threshold``
     - 2.0
     - Inlier threshold (raw, not squared)
   * - ``confidence``
     - 0.99
     - RANSAC confidence
   * - ``max_iterations``
     - 5 000
     - Maximum iterations
   * - ``spatial_coherence_weight``
     - 0.975
     - Weight for graph-cut spatial coherence (0–1)
   * - ``cell_number``
     - 16
     - Grid cells per dimension for neighbourhood graph
   * - ``do_local_optimization``
     - ``True``
     - Enable GC-RANSAC local optimisation
   * - ``do_graph_cut``
     - ``True``
     - Enable graph-cut inlier selection
   * - ``do_final_least_squares``
     - ``True``
     - Final refinement on inliers

**When to use:** When spatially coherent inlier selection matters — e.g.
natural images with clustered outlier distributions.


MAGSAC++
--------

**MAGSAC++** [Barath2020]_ — threshold-free scoring via sigma consensus.
Instead of a hard inlier/outlier boundary, MAGSAC++ marginalises over
a range of possible noise levels, making it less sensitive to the choice
of threshold.

.. code-block:: python

   solver = homlib.Wadenback3DV2026OneSided()
   mag = homlib.MAGSAC(solver, sigma_max=10.0)
   estimate, stats = mag.estimate(x, y)

Requires building with ``-DHOMLIB_USE_MAGSAC=ON``. Check availability:

.. code-block:: python

   if homlib._HAS_MAGSAC:
       mag = homlib.MAGSAC(solver, sigma_max=10.0)

**Key settings** (via :class:`~homlib.MAGSACSettings`):

.. list-table::
   :header-rows: 1
   :widths: 35 15 50

   * - Parameter
     - Default
     - Description
   * - ``sigma_max``
     - 10.0
     - Maximum noise level for sigma consensus
   * - ``confidence``
     - 0.99
     - RANSAC confidence
   * - ``max_iterations``
     - 5 000
     - Maximum iterations
   * - ``partition_number``
     - 5
     - Partitions in the sigma consensus integral

**When to use:** When you want to avoid tuning a precise inlier threshold.
The sigma consensus scoring is more forgiving of threshold mis-specification.


SupeRANSAC
----------

**SupeRANSAC** [Barath2025]_ — a highly modular backend where the sampler,
scorer, and local optimiser are independently configurable.

.. code-block:: python

   import homlib

   solver = homlib.Wadenback3DV2026OneSided()

   sa = homlib.SupeRANSAC(
       solver,
       inlier_threshold=0.07,
       scoring_type=homlib.SuperansacScoringType.MAGSAC,
       sampler_type=homlib.SuperansacSamplerType.PROSAC,
       lo_type=homlib.SuperansacLocalOptType.IRLS,
   )
   estimate, stats = sa.estimate(x, y)

Requires building with ``-DHOMLIB_USE_SUPERANSAC=ON``. Check availability:

.. code-block:: python

   if homlib._HAS_SUPERANSAC:
       sa = homlib.SupeRANSAC(solver, inlier_threshold=0.07)

**Configurable components:**

- **Scoring:** ``SuperansacScoringType`` — RANSAC, MSAC, MAGSAC, MINPRAN,
  ACRANSAC, and others.
- **Sampling:** ``SuperansacSamplerType`` — Uniform, PROSAC, NAPSAC,
  ProgressiveNAPSAC, ImportanceSampler, ARSampler.
- **Local optimisation:** ``SuperansacLocalOptType`` — None, LSQ, IRLS,
  NestedRANSAC, GCRANSAC, IteratedLMEDS, CrossValidation.

**When to use:** Research and benchmarking, or when you need a specific
combination of RANSAC components not available in the other backends.


.. _threshold_conventions:

Threshold conventions
---------------------

The backends use different threshold conventions:

- **LOMSAC**: ``inlier_threshold`` (raw, same as other backends) or
  ``squared_inlier_threshold`` (squared).
  ``inlier_threshold = pixel_threshold / focal_length``
- **GC-RANSAC**: raw threshold in normalised coordinates.
  ``threshold = pixel_threshold / focal_length``
- **MAGSAC++**: sigma maximum (upper bound on noise sigma).
  ``sigma_max`` in normalised coordinates.
- **SupeRANSAC**: raw threshold in normalised coordinates.
  ``inlier_threshold = pixel_threshold / focal_length``

All coordinates must be centred and normalised before calling ``estimate()``.
See :ref:`solver_guide` for details on coordinate conventions.


Choosing a backend
------------------

- **Default / production:** :class:`~homlib.LOMSAC` — always available,
  deterministic with ``random_seed``, well-tested.
- **Threshold-free scoring:** :class:`~homlib.MAGSAC` — when you do not want
  to fine-tune an inlier threshold.
- **Spatial coherence:** :class:`~homlib.GCRANSAC` — when outlier distributions
  are spatially structured.
- **Research / customisation:** :class:`~homlib.SupeRANSAC` — maximum
  flexibility for ablation studies.

.. seealso::

   :ref:`api` for the full API reference, :ref:`references` for the original
   papers, and the :ref:`sphx_glr_gallery` for worked examples.
