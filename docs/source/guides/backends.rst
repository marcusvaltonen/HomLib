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

**MAGSAC++** [Barath2019]_ [Barath2020]_ — threshold-free scoring via sigma consensus.
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
       scoring=homlib.SuperansacScoringType.MAGSAC,
       sampler=homlib.SuperansacSamplerType.PROSAC,
       local_optimization=homlib.SuperansacLocalOptType.IRLS,
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


.. _confidence_scores:

Confidence scores and sampler selection
---------------------------------------

When correspondences have associated confidence scores (e.g. from a learned
matcher like LoMa, SuperGlue, or LightGlue), they can be passed to the RANSAC
backends to guide sampling towards more reliable matches. This is supported by
GC-RANSAC, MAGSAC, and SupeRANSAC.

.. code-block:: python

   # scores: 1-D array of per-match confidences in [0, 1], length N,
   # sorted in descending order (highest confidence first).
   gc = homlib.GCRANSAC(solver, threshold=3.0,
                        sampler=homlib.GCSamplerType.PROSAC)
   estimate, stats = gc.estimate(x, y, scores=scores)

**Available samplers for GC-RANSAC and MAGSAC** (``homlib.GCSamplerType``):

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Sampler
     - Description
   * - ``Uniform``
     - Standard uniform random sampling (default, ignores scores).
   * - ``PROSAC``
     - Progressively expands the sampling pool from the highest-scored
       matches outward. Requires input sorted by confidence but does not
       use the score values directly.
   * - ``ARSampler``
     - Adaptive reordering sampler. Uses the actual score values to weight
       the sampling distribution.
   * - ``ImportanceSampler``
     - Samples proportionally to the score values.

**Available samplers for SupeRANSAC** (``homlib.SuperansacSamplerType``):
the same set as above, plus ``NAPSAC`` and ``ProgressiveNAPSAC``
(neighbourhood-based samplers).

Scores are **optional**: omitting them uses Uniform sampling. If a sampler
that needs scores is selected without providing them, a ``RuntimeError`` is
raised (GC-RANSAC, MAGSAC) or the sampler falls back to uniform weights
(SupeRANSAC ImportanceSampler).

**Score normalization.** Scores should be normalized to [0, 1]. Two common
approaches:

- **Real scores** — min-max normalization of the raw matcher confidences.
  Preserves the relative magnitudes between matches, which benefits
  samplers that use the actual values (ARSampler, ImportanceSampler).
- **Linear scores** — replace the values with
  ``np.linspace(1.0, 0.01, N)``, preserving the sorted order but
  discarding the actual values.  This is the approach used by Barath
  et al. in the SupeRANSAC and GC-RANSAC codebases: because PROSAC
  only uses the rank order (not the magnitudes), the linear ramp is
  sufficient and avoids sensitivity to the matcher's score distribution.

.. code-block:: python

   import numpy as np

   # Real scores (from a learned matcher, already sorted descending)
   scores_raw = pair["scores"]
   s_min, s_max = scores_raw.min(), scores_raw.max()
   scores_real = (scores_raw - s_min) / (s_max - s_min)

   # Linear scores (Barath et al. convention)
   # Works well with PROSAC, which only cares about ordering.
   N = x.shape[1]
   scores_linear = np.linspace(1.0, 0.01, N)


**Sampler trade-offs.**

- **Uniform** is the fastest sampler (no overhead), but ignores match
  quality entirely — convergence depends on the inlier ratio.
- **PROSAC** adds minimal overhead and converges faster than Uniform
  when good matches are at the front of the sorted list, because it
  draws early samples from the top-ranked correspondences.  It is
  generally the best default when scores are available.
- **ARSampler** and **ImportanceSampler** use the actual score values
  to weight the sampling distribution.  They can outperform PROSAC
  when the score magnitudes are informative (e.g. learned matcher
  confidences), but add more computational overhead per iteration.

On real-world data (e.g. HPatches), PROSAC with real scores typically
matches or slightly improves Uniform accuracy while needing fewer
iterations to converge.  The wall-clock difference is small because
the per-iteration cost dominates, but PROSAC is particularly helpful
at low inlier ratios where Uniform needs many more iterations to find
a good sample.

See the :ref:`sampler_scores_example` for a runnable example.


.. _backend_settings:

Backend settings
----------------

Each backend exposes its full set of configurable parameters as keyword
arguments. Settings are also accessible via namespace aliases for clarity:

.. code-block:: python

   # GC-RANSAC
   gc = homlib.GCRANSAC(solver,
       threshold=3.0,
       confidence=0.99,
       max_iterations=5000,
       min_iterations=50,
       do_local_optimization=True,
       sampler=homlib.GCRANSAC.SamplerType.PROSAC,
       non_minimal_strategy=homlib.GCRANSAC.NonMinimalStrategy.SelectAndRefine,
   )

   # MAGSAC (shares GCSamplerType with GC-RANSAC)
   mg = homlib.MAGSAC(solver,
       sigma_max=3.0,
       sampler=homlib.MAGSAC.SamplerType.ARSampler,
   )

   # SupeRANSAC — all sub-settings are exposed
   sa = homlib.SupeRANSAC(solver,
       inlier_threshold=3.0,
       sampler=homlib.SupeRANSAC.SamplerType.ARSampler,
       scoring=homlib.SuperansacScoringType.MSAC,
       local_optimization=homlib.SuperansacLocalOptType.GCRANSAC,
   )

**Namespace aliases.** Each backend class provides aliases to its associated
types, so you can write ``homlib.GCRANSAC.SamplerType`` instead of
``homlib.GCSamplerType``:

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Backend
     - Alias
     - Resolves to
   * - GC-RANSAC
     - ``GCRANSAC.Settings``
     - ``GCRANSACSettings``
   * -
     - ``GCRANSAC.SamplerType``
     - ``GCSamplerType``
   * -
     - ``GCRANSAC.NonMinimalStrategy``
     - ``NonMinimalStrategy``
   * - MAGSAC
     - ``MAGSAC.Settings``
     - ``MAGSACSettings``
   * -
     - ``MAGSAC.SamplerType``
     - ``GCSamplerType``
   * -
     - ``MAGSAC.NonMinimalStrategy``
     - ``NonMinimalStrategy``
   * - SupeRANSAC
     - ``SupeRANSAC.Settings``
     - ``SuperansacSettings``
   * -
     - ``SupeRANSAC.SamplerType``
     - ``SuperansacSamplerType``
   * -
     - ``SupeRANSAC.ScoringType``
     - ``SuperansacScoringType``
   * -
     - ``SupeRANSAC.LocalOptType``
     - ``SuperansacLocalOptType``

**SupeRANSAC sub-settings.** SupeRANSAC exposes nested settings structs for
fine-grained control:

.. code-block:: python

   settings = homlib.SuperansacSettings()
   settings.ar_sampler_settings.estimator_variance = 0.5
   settings.lo_settings.spatial_coherence_weight = 0.14
   settings.lo_settings.graph_cut_number = 3
   settings.degensac_settings.enabled = True
   settings.neighborhood_settings.neighborhood_size = 10.0

   sa = homlib.SupeRANSAC(solver, settings=settings)

The sub-settings structs are:

- ``SuperansacARSamplerSettings`` — AR sampler variance and randomness.
- ``SuperansacLOSettings`` — local optimisation iterations, graph-cut
  settings, spatial coherence weight.
- ``SuperansacDegensacSettings`` — degenerate model handling.
- ``SuperansacNeighborhoodSettings`` — neighbourhood graph parameters.

**Non-minimal strategy** (GC-RANSAC and MAGSAC). Controls how the
non-minimal sample is used for model refinement:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Strategy
     - Description
   * - ``Lomsac``
     - LO-MSAC: iterative refinement on inlier set.
   * - ``SelectAndRefine``
     - Select the best model, then refine (default).
   * - ``MultiSample``
     - Evaluate multiple non-minimal samples.


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


.. _backend_availability:

Platform availability
---------------------

.. list-table::
   :header-rows: 1
   :widths: 20 20 20 20 20

   * - Backend
     - Linux
     - macOS
     - Windows
     - PyPI wheels
   * - LOMSAC
     - Yes
     - Yes
     - Yes
     - All platforms
   * - GC-RANSAC
     - Yes
     - Yes
     - Yes
     - —
   * - MAGSAC++
     - Yes
     - Yes
     - Yes
     - —
   * - SupeRANSAC
     - Yes
     - Yes
     - Yes
     - —

**Pre-built wheels** (``pip install homlib``) include **LOMSAC only**.

**Building from source** enables all four backends when the dependencies are
installed:

.. code-block:: bash

   # Linux (Debian/Ubuntu)
   sudo apt-get install libboost-dev libboost-system-dev libopencv-dev libgflags-dev libomp-dev

   # macOS
   brew install boost opencv gflags libomp

   # Windows (vcpkg)
   vcpkg install boost opencv gflags openmp --triplet x64-windows

   # Build with all backends
   pip install .

To check which backends are available at runtime:

.. code-block:: python

   import homlib
   print("SupeRANSAC:", homlib._HAS_SUPERANSAC)
   print("GC-RANSAC: ", homlib._HAS_GCRANSAC)
   print("MAGSAC++:  ", homlib._HAS_MAGSAC)


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
   papers, and the :ref:`examples` for worked examples.
