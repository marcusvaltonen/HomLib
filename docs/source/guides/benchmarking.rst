.. _benchmarking_guide:

Benchmarking
============

HomLib's benchmarks compare its solvers and robust estimators on real image
pairs (HPatches) and on synthetic problems. Use them to choose a solver for
your data, and to measure the effect of a change in HomLib: run a benchmark
before and after the change, and compare the two runs.

The scripts are in ``python/experiments/``. They need the ``homlib`` Python
package (see :ref:`installation_users`), and a few packages of their own:

.. code-block:: bash

   pip install opencv-python tqdm

The commands below run from the ``python/`` directory of the repository.


HPatches benchmark
------------------

What it measures
^^^^^^^^^^^^^^^^

`HPatches <https://github.com/hpatches/hpatches-dataset>`__ (Balntas et al.,
CVPR 2017) has 116 sequences of 6 images of planar scenes. In the 59
*viewpoint* sequences, the camera moves; in the 57 *illumination* sequences,
the camera stays still and the lighting changes. Each sequence gives the
ground-truth homography from its first image to each of the 5 others: 580
pairs in all.

The HPatches images have no lens distortion. The benchmark therefore measures
what each solver gets right on real matches, with their noise and outliers,
and what the distortion solvers lose by estimating parameters that are zero.

The benchmark runs a robust estimator (LOMSAC by default) with every point
solver on every pair, 3 times per pair, and reports three metrics,
separately for the viewpoint and the illumination pairs:

**Inliers found**
   The share of the true inliers that the estimator returns as inliers. The
   *true inliers* are the matches within the inlier threshold (3 px by
   default) of the ground-truth homography. This metric ignores false
   inliers: the corner error accounts for them. Pairs with no true inlier are
   left out, since no solver can find any.

**AUC@1, AUC@3, AUC@5, AUC@10**
   The area under the cumulative curve of the corner error, up to 1, 3, 5 and
   10 px, in percent. The *corner error* of a pair maps the four corners of the
   first image with the estimated model and with the ground-truth homography,
   and averages the distances between the two. The estimated model is the full
   one: the corners are undistorted with the estimated distortion of the first
   image, mapped with the estimated homography, and distorted with the
   estimated distortion of the second image. A pair with an error of 0 scores
   100; a pair with an error above the threshold scores 0. Wadenbäck et al.
   [Wadenback2026]_ report the AUC at 1, 3 and 5 px on HPatches, as "mAA
   (AUC)".

**Time**
   The median time of a run of the robust estimator, in milliseconds,
   measured from Python.

Step 1: compute the matches
^^^^^^^^^^^^^^^^^^^^^^^^^^^

``hpatches_features.py`` downloads HPatches if needed (a large archive),
matches the 580 pairs, and saves the matched points with the ground-truth
homographies to a ``.npz`` file:

.. code-block:: bash

   python experiments/hpatches_features.py --data_dir data --output data/hpatches_sift.npz

Matching is a separate step, so that every later run compares the solvers on
the same matches. The script has three matchers:

``--features sift`` (default)
   OpenCV SIFT, matched with Lowe's ratio test. It needs no learned model, and
   is easy to reproduce. With OpenCV 4.8 or later, it uses SIFT's precise
   upscaling, which the accuracy of the keypoints needs.

``--matcher lightglue``, with ``--features superpoint``, ``disk``, ``aliked`` or ``sift``
   These features, matched with
   `LightGlue <https://github.com/cvg/LightGlue>`__ (Lindenberger et al., ICCV
   2023). For ``sift``, the script gives LightGlue the same OpenCV SIFT with
   precise upscaling: LightGlue's own SIFT lacks it, which shifts its points by
   about 0.1 px.

``--features loma``
   `LoMa <https://github.com/davnords/LoMa>`__ (Nordström et al., ECCV 2026),
   which detects, describes and matches the points itself.

LightGlue and LoMa install from their GitHub repositories (``pip install -e .``
in a clone; the PyPI package named ``loma`` is unrelated), and run much faster
on a GPU. SuperPoint's weights are for non-commercial use only.

Wadenbäck et al. [Wadenback2026]_ evaluate on HPatches with SuperPoint and
LightGlue, 1024 keypoints per image, and the shorter side of the images
resized to 480 px. These settings give the same matches:

.. code-block:: bash

   python experiments/hpatches_features.py --data_dir data --features superpoint \
       --matcher lightglue --num_features 1024 --resize 480 \
       --output data/hpatches_superpoint_lightglue.npz

The saved points are in pixels of the original images, centred on each
image's centre, which is HomLib's distortion centre; the homography is in these
coordinates. OpenCV and LightGlue put the origin of the pixel coordinates at
the centre of the top-left pixel, and LoMa at the top-left corner of the image:
the script shifts LoMa's points by -0.5 px, to OpenCV's convention. Since the
points are not normalised, the inlier thresholds of the benchmark are in pixels.

Step 2: run the solvers
^^^^^^^^^^^^^^^^^^^^^^^

``hpatches_benchmark.py`` runs the solvers on a matches file and prints a table
of the three metrics:

.. code-block:: bash

   python experiments/hpatches_benchmark.py data/hpatches_sift.npz --output sift.json

Its options:

``--solvers NAME ...``
   The solvers to run, by their :class:`~homlib.SolverId` names (default:
   every point solver).

``--backend NAME``
   The robust estimator: ``lomsac`` (default), ``gcransac``, ``magsac`` or
   ``superansac``. GC-RANSAC, MAGSAC++ and SupeRANSAC need HomLib built with
   them (see :ref:`backend_availability`). MAGSAC++ has no inlier threshold:
   it gets ``--threshold`` as its largest noise level, ``sigma_max``.

``--threshold PX``
   The inlier threshold, in pixels (default: 3).

``--seeds N``
   The number of runs per pair (default: 3). LOMSAC gets a different random
   seed for each run, so that its runs are reproducible; the other estimators
   draw their seeds themselves.

``--iterations N``
   Runs the robust estimator for exactly N iterations, instead of stopping
   when it is confident. Every solver then gets the same budget of
   iterations. Several runs with different budgets give the accuracy as a
   function of the time, as in Fig. 6 of Wadenbäck et al. [Wadenback2026]_.

``--output FILE``
   Saves the results of every pair and seed to a JSON file.

``--compare A B``
   Prints the differences B - A between two saved runs, and exits. The
   differences are paired: for each solver and each pair present in both
   runs, the script takes the difference of the share of the true inliers
   found and of AUC@5, averaged over the seeds. It reports their mean, with
   its standard error over the pairs, and the ratio of the median times.
   A difference smaller than about twice its standard error is not
   significant.

To measure the effect of a change in HomLib:

.. code-block:: bash

   python experiments/hpatches_benchmark.py data/hpatches_sift.npz --output before.json
   # ... change HomLib, and reinstall it ...
   python experiments/hpatches_benchmark.py data/hpatches_sift.npz --output after.json
   python experiments/hpatches_benchmark.py --compare before.json after.json

Two other scripts evaluate on HPatches: ``hpatches.py`` with
`glue-factory <https://github.com/cvg/glue-factory>`__, and
``hpatches_loma.py`` with LoMa matches computed on the fly.


Synthetic benchmark
-------------------

What it measures
^^^^^^^^^^^^^^^^

``synthetic_benchmark.py`` solves problems from
:func:`~homlib.generate_problem_instance`. Each problem is a random plane seen
by two cameras with a random relative pose and a 70° field of view. The
generator distorts the points of each image with the division model, with a
random parameter between -0.2 and -0.01, on the side or sides that its
:class:`~homlib.DistortionCase` sets, and returns the points, the
ground-truth homography and the distortion parameters. The benchmark adds
Gaussian noise to the points, given in pixels for a focal length of 1000 px.

Unlike HPatches, these problems have distortion and an exact solution. The
benchmark measures each solver on the distortion cases it models (``NO_DISTORTION``
and the case it is made for, or all five for the general two-sided solvers),
with the experiments of Wadenbäck et al. [Wadenback2026]_:

**Numerical stability** (their Fig. 3)
   Without noise, how often the solver misses the true solution: none of its
   solutions is within 1e-6 of it, in homography error and in the sum of the
   errors of the two distortion parameters. A stable solver misses almost
   never.

**Sensitivity to noise** (their Fig. 4)
   With noise of 0.5, 1 and 2 px, the median error of the solver's best
   solution, over the problems where it returns a solution.

**No solution**
   How often the solver returns no solution at all.

**Time**
   The median time of a call to the solver, in microseconds, measured from
   Python: it includes the overhead of the bindings.

Each solver gets minimal samples: 4 points for the solvers without distortion,
5 or 6 for the others. By default, the benchmark solves 2000 problems per
solver, case and noise level:

.. code-block:: bash

   python experiments/synthetic_benchmark.py --output synthetic.json

``--problems N``, ``--noise PX ...`` and ``--solvers NAME ...`` change these
settings.

LOMSAC with outliers
^^^^^^^^^^^^^^^^^^^^

``--outliers`` runs another experiment, as in Fig. 5 of Wadenbäck et al.
[Wadenback2026]_: LOMSAC on problems with 200 points, 0.5 px of noise and the
given shares of outliers, with a 5 px inlier threshold. It reports the share of
the true inliers that LOMSAC finds, and the median time of a run. Each solver
solves problems of the distortion case it is made for, so the rows of solvers
made for different cases don't compare.

.. code-block:: bash

   python experiments/synthetic_benchmark.py --outliers 0.2 0.4 0.6 0.8 --output outliers.json

Without ``--iterations``, LOMSAC stops when it is confident, or after 10000
iterations. ``--iterations N`` runs exactly N iterations, so that every solver
gets the same budget. By default, this experiment solves 200 problems per
solver and share of outliers.

Comparing two runs
^^^^^^^^^^^^^^^^^^

``--output`` and ``--compare A B`` work as in the HPatches benchmark. Here the
two runs don't solve the same problems, since
:func:`~homlib.generate_problem_instance` seeds itself with the time: the
differences are statistical, and the script reports them with their standard
errors.


Regression test
---------------

``python/tests/test_hpatches.py`` runs every point solver in LOMSAC on 8 HPatches
pairs (4 viewpoint, 4 illumination), with OpenCV SIFT matches, a 3 px inlier
threshold and 3 seeds. It checks that each solver finds, on average, at least
95% of the true inliers. The two-sided solvers have lower bounds, from 75% to
85%: they are weaker on the illumination pairs, whose homography is the
identity, which any equal distortion parameters fit. The test runs in a few
seconds, with the other tests:

.. code-block:: bash

   python -m pytest tests/test_hpatches.py

Test data
^^^^^^^^^

The matches of the 8 pairs are not in the repository. They are the asset
``hpatches_sift_subset.npz`` of the GitHub release ``test-data-v1``, published
next to ``hpatches_sift_subset.npz.sha256``, its SHA-256 checksum. The test
looks for the file, in this order:

1. ``$HOMLIB_TEST_DATA/hpatches_sift_subset.npz``, if the environment variable
   ``HOMLIB_TEST_DATA`` names a directory that contains it;
2. a cached copy in ``~/.cache/homlib/`` (or ``$XDG_CACHE_HOME/homlib/``);
3. a download from the release, checked against the SHA-256 file, then
   cached. ``HOMLIB_TEST_DATA_URL`` replaces the URL of the release asset, for
   example to test a fork.

Without the file, the test is skipped. If ``HOMLIB_REQUIRE_TEST_DATA`` is
set, as in CI, the test fails instead, so that a missing file doesn't go
unnoticed.

Publishing the test data (maintainers)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The release must exist before CI runs the test. To write the file from the
full SIFT matches and publish it with its checksum:

.. code-block:: bash

   python experiments/hpatches_features.py --matches data/hpatches_sift.npz \
       --export_subset hpatches_sift_subset.npz
   shasum -a 256 hpatches_sift_subset.npz > hpatches_sift_subset.npz.sha256
   gh release create test-data-v1 hpatches_sift_subset.npz hpatches_sift_subset.npz.sha256 \
       -R marcusvaltonen/HomLib --latest=false \
       --title "Test data" --notes "Data files for the tests (not a HomLib release)"

To change the data, publish it under a new tag (``test-data-v2``) and update
``SUBSET_URL`` in the test: existing checkouts and caches keep the old file.
HPatches is under the MIT license, which allows publishing these matches.


Results
-------

The table gives the HPatches results of five solvers with OpenCV SIFT matches,
LOMSAC, a 3 px threshold and 3 seeds per pair, on an AWS g4dn.xlarge host
(Ubuntu 24.04). The times depend on the host: compare them within a run.

.. list-table::
   :header-rows: 1
   :widths: 34 16 14 9 9 10

   * - Solver
     - Sequences
     - Inliers found
     - AUC@3
     - AUC@5
     - Time (ms)
   * - :class:`~homlib.DLT`
     - viewpoint
     - 93.35%
     - 40.68
     - 52.90
     - 8.50
   * -
     - illumination
     - 96.30%
     - 64.51
     - 75.58
     - 6.46
   * - :class:`~homlib.GuoPRL2019`
     - viewpoint
     - 93.29%
     - 40.65
     - 52.95
     - 5.74
   * -
     - illumination
     - 96.22%
     - 64.52
     - 75.58
     - 3.72
   * - :class:`~homlib.FitzgibbonCVPR2001OneSided`
     - viewpoint
     - 93.03%
     - 20.53
     - 31.93
     - 19.24
   * -
     - illumination
     - 95.62%
     - 60.13
     - 72.52
     - 13.47
   * - :class:`~homlib.Wadenback3DV2026OneSided`
     - viewpoint
     - 92.71%
     - 20.88
     - 32.47
     - 13.47
   * -
     - illumination
     - 95.56%
     - 59.87
     - 72.08
     - 8.61
   * - :class:`~homlib.Wadenback3DV2026TwoSided`
     - viewpoint
     - 92.94%
     - 14.55
     - 25.23
     - 18.14
   * -
     - illumination
     - 91.66%
     - 19.30
     - 29.70
     - 8.39

All five solvers find most of the true inliers. The corner error separates
them: on these undistorted images, the solvers without distortion are the most
accurate, and the more distortion parameters a solver estimates, the lower
its AUC. The two-sided solver loses most on the illumination pairs, where any
equal distortion parameters fit the identity homography.

The full results are in the description of
`pull request #44 <https://github.com/marcusvaltonen/HomLib/pull/44>`__: every
point solver on HPatches with SIFT, SuperPoint + LightGlue and LoMa matches,
the synthetic experiments of Figs. 3, 4 and 5, and the AUC@5 on HPatches for a
fixed number of iterations (Fig. 6).
