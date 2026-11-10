Synthetic benchmark
===================

``synthetic_benchmark.py`` solves random problems from
:func:`~homlib.generate_problem_instance`: a plane seen by two cameras with a
70° field of view, with division-model distortion (a parameter between -0.2
and -0.01) on the side or sides of each :class:`~homlib.DistortionCase`, and
Gaussian noise given in pixels for a focal length of 1000 px. Unlike HPatches,
these problems have distortion and an exact solution.

Each solver gets minimal samples and solves the distortion cases it models.
The *error* of a solution is the larger of the homography error (the Frobenius
norm of the difference with the true homography, both scaled to a norm of 1)
and of the sum of the errors of the distortion parameters; the *best solution*
is the one with the smallest error. The benchmark reports, as the experiments
of Wadenbäck et al. [Wadenback2026]_:

**Missed** (their Fig. 3)
   Without noise, the share of the problems where the best solution has an
   error above 1e-6: a stable solver almost never misses.

**Median error** (their Fig. 4)
   With 0.5, 1 and 2 px of noise, the median error of the best solution, over the
   problems with a solution.

**No solution**
   The share of the problems where the solver returns none.

**Time**
   The median time of a call to the solver, from Python.

Running it
----------

.. code-block:: bash

   python benchmark/synthetic/synthetic_benchmark.py --output synthetic.json

By default it solves 2000 problems per solver, case and noise level
(``--problems``, ``--noise`` and ``--solvers`` change them). The problems are
seeded (``--seed``), so two runs solve the same problems and give the same
results on the same platform.

``--outliers 0.2 0.4 0.6 0.8`` runs another experiment, as Fig. 5 of Wadenbäck
et al.: LOMSAC on problems with 200 points, 0.5 px of noise and these shares
of outliers, with a 5 px threshold; it reports the share of the true inliers
found and the median time of a run. ``--problems`` then defaults to 200,
``--noise`` doesn't apply, and ``--iterations N`` and ``--backend`` (see the
:ref:`common options <benchmark_common_options>`) are available.

.. _synthetic_results:

Results
-------

Five solvers with the default settings, each on the distortion case it is made
for, on an AWS g4dn.xlarge host (Linux). *Missed* is for the noise-free
problems, *No solution* for 1 px of noise, and the median errors for 0.5, 1
and 2 px; the times include the Python bindings.

.. list-table::
   :header-rows: 1

   * - Solver
     - Case
     - Missed
     - No solution (1 px)
     - Median error (0.5, 1, 2 px)
     - Time
   * - :class:`~homlib.DLT`
     - ``NO_DISTORTION``
     - 0.00%
     - 1.30%
     - 0.038, 0.078, 0.154
     - 25 µs
   * - :class:`~homlib.GuoPRL2019`
     - ``NO_DISTORTION``
     - 0.00%
     - 1.30%
     - 0.038, 0.078, 0.154
     - 7 µs
   * - :class:`~homlib.FitzgibbonCVPR2001OneSided`
     - ``ONE_SIDED_LEFT``
     - 0.00%
     - 0.00%
     - 0.076, 0.149, 0.345
     - 33 µs
   * - :class:`~homlib.Wadenback3DV2026OneSided`
     - ``ONE_SIDED_LEFT``
     - 0.00%
     - 3.45%
     - 0.065, 0.132, 0.270
     - 9 µs
   * - :class:`~homlib.Wadenback3DV2026TwoSided`
     - ``TWO_SIDED``
     - 0.30%
     - 1.60%
     - 0.139, 0.265, 0.507
     - 21 µs

DLT and Guo give the same homography from the same 4 points (it is unique), so
the same errors; Guo's closed form is faster. The more parameters a solver
estimates, the larger its error for the same noise. With noise, some solvers
sometimes return no solution: for the 4-point solvers, when their sample check
rejects a sample with three nearly collinear points.
