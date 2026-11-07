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

The benchmark runs a robust estimator (:class:`~homlib.LOMSAC` by default,
HomLib's LO-RANSAC: see :ref:`backend_guide`) with every *point solver* on
every pair, 3 times per pair. The point solvers are the solvers that take point
matches alone, as opposed to the solvers of affine or orientation features (see
:ref:`feature_types_overview`). The benchmark reports three metrics, separately
for the viewpoint and the illumination pairs:

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
   estimated distortion of the second image (division model, see
   :ref:`solver_guide`). Each run with a corner error :math:`e` scores
   :math:`\max(0, 1 - e / t)` at the threshold :math:`t`, and the AUC@t is the
   mean of these scores over the pairs and the seeds, in percent: an error of 0
   scores 100, an error above the threshold scores 0. A failed estimate,
   whose model doesn't map the corners to finite points, counts as an infinite
   error, and scores 0. Wadenbäck et al. [Wadenback2026]_ report the AUC at 1,
   3 and 5 px on HPatches, as "mAA (AUC)".

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
   upscaling, which the accuracy of the keypoints needs. By default, it detects
   4000 keypoints per image (``--num_features``) and keeps the matches with a
   ratio below 0.8 (``--ratio``): the `Results`_ below use these defaults.

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
resized to 480 px. These settings follow their setup:

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

With ``--resize``, the features are detected in the resized images, but the
points stay in pixels of the original images, and so do the inlier thresholds
and the corner errors of the benchmark by default. glue-factory's HPatches
evaluation measures them in pixels of the resized images, and so, probably, do
Wadenbäck et al.: ``--protocol gluefactory`` (see `glue-factory's protocol`_)
measures them as it does, so that the AUCs compare with theirs.

The matches files also store the sizes of the two images of each pair: the
images of a viewpoint sequence don't all have the same size. Older matches
files, such as the regression test's subset, lack the second image's size;
``--protocol gluefactory`` then reads it from the HPatches images
(``--data_dir``).

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
   The inlier threshold, in pixels (default: 3): of the original images, or of
   the resized images with ``--protocol gluefactory``.

``--protocol original|gluefactory``
   ``original`` (default) measures in pixels of the original images, on the
   580 pairs. ``gluefactory`` follows glue-factory's HPatches evaluation (see
   `glue-factory's protocol`_). With matches files that lack the sizes of the
   second images, it needs ``--data_dir``, the directory that contains
   ``hpatches-sequences-release``, as for ``hpatches_features.py``.

``--seeds N``
   The number of runs per pair (default: 3). LOMSAC gets a different random
   seed for each run, so that its runs are reproducible on the same platform;
   the other estimators draw their seeds themselves. Across platforms, the
   results differ slightly, because the standard library's
   ``std::uniform_int_distribution`` is implementation-defined: libc++ (macOS)
   and libstdc++ (Linux) draw different samples from the same seed. For
   example, ``DLT`` finds 93.35% of the true inliers of the viewpoint pairs
   of the SIFT matches on Linux, and 93.14% on macOS.

``--iterations N``
   Runs the robust estimator for exactly N iterations, instead of stopping
   when it is confident. Every solver then gets the same budget of
   iterations. Several runs with different budgets give the accuracy as a
   function of the time, as in Fig. 6 of Wadenbäck et al. [Wadenback2026]_.

``--output FILE``
   Saves the settings (the name of the matches file, ``--threshold``,
   ``--backend``, ``--seeds``, ``--iterations`` and ``--protocol``) and the
   results of every pair and seed to a JSON file. Files without
   ``--protocol`` in their settings count as ``original``.

``--compare A B``
   Prints the differences B - A between two saved runs, and exits. The
   differences are paired: for each solver and each pair present in both
   runs, the script takes the difference of the share of the true inliers
   found and of AUC@5, averaged over the seeds. It reports their mean, with
   its standard error over the pairs, and the ratio of the median times,
   which has no standard error. A difference smaller than about twice its
   standard error is not significant. The script warns if the two runs have
   different settings, or if a file has none (files saved by an older
   version of the script). For two runs with ``--protocol gluefactory``, it
   also compares all the pairs together.

To measure the effect of a change in HomLib:

.. code-block:: bash

   python experiments/hpatches_benchmark.py data/hpatches_sift.npz --output before.json
   # ... change HomLib, and reinstall it ...
   python experiments/hpatches_benchmark.py data/hpatches_sift.npz --output after.json
   python experiments/hpatches_benchmark.py --compare before.json after.json

Another script, ``hpatches_loma.py``, evaluates on HPatches with LoMa matches
computed on the fly.

glue-factory's protocol
^^^^^^^^^^^^^^^^^^^^^^^

`glue-factory <https://github.com/cvg/glue-factory>`__'s HPatches evaluation
(``gluefactory/eval/hpatches.py``, with the dataset of
``gluefactory/datasets/hpatches.py``) is the protocol of the HPatches results
in glue-factory's README. ``--protocol gluefactory`` follows it:

- **Resized images.** glue-factory resizes each image so that its shorter side
  is 480 px, the other side rounded down: a 1000 × 750 image becomes
  640 × 480. The two images of a pair are resized separately. The benchmark
  scales the points, the ground-truth homography and the image corners
  accordingly, so that the inlier threshold, the true inliers and the corner
  error are in pixels of the resized images. The features aren't detected
  again: the matches are those of ``hpatches_features.py``, detected in the
  original images (or in resized ones, with its ``--resize 480``), and only
  rescaled. Resizing can also enlarge an image: the first images of 65 of the
  580 pairs have a shorter side under 480 px, and the scale factors range
  from 0.23 to 1.26.
- **Ground truth.** glue-factory applies the homographies of HPatches to
  coordinates with their origin at the top-left corner of the image, whereas
  ``hpatches_features.py`` applies them with the origin at the centre of the
  top-left pixel, as OpenCV does. Which convention HPatches itself uses is
  undetermined. ``--protocol gluefactory`` follows glue-factory's: the two
  ground truths map the corners 0.08 px apart on the median pair, but more
  than 1 px apart on 1.7% of the pairs, and up to 5 px (``v_charing/3``).
- **540 pairs.** glue-factory leaves out 8 sequences of large images
  (``i_contruction``, ``i_crownnight``, ``i_dc``, ``i_pencils``,
  ``i_whitebuilding``, ``v_artisans``, ``v_astronautis``, ``v_talent``), which
  leaves 52 illumination and 56 viewpoint sequences.
- **Every pair counts.** The pairs with no true inlier count in the AUC, but
  not in the share of the true inliers found, which is undefined for them. An
  estimate with no inlier is a failure, with an infinite corner error: this
  stands in for glue-factory's failure, when its estimator returns no model.
- **Corner error and AUC.** Its corner error is the benchmark's: the mean
  distance between the four outer corners of the first image mapped with the
  estimated and with the ground-truth homography (the benchmark maps them with
  the full estimated model, distortion included). It reports the AUC at 1, 3
  and 5 px over all the pairs: the table gets a row for all the pairs, besides
  the illumination and the viewpoint pairs. It integrates the AUC with the
  trapezoidal rule, which adds less than :math:`50 / n` points to the
  benchmark's AUC, for :math:`n` runs: less than 0.1 points on 540 pairs.
- **Robust estimator and threshold.** The robust estimator stays HomLib's
  (LOMSAC by default, or ``--backend``), whereas glue-factory runs a
  homography RANSAC of its own: OpenCV's ``findHomography`` with a 0.5 px
  threshold in its official configurations of LightGlue, or PoseLib with a
  1 px threshold in the pipeline's default, without a configuration.
  ``--threshold`` sets the threshold, in pixels of the resized images. There
  is no threshold sweep: glue-factory's ``eval.ransac_th=-1`` runs 0.5, 1,
  1.5, 2, 2.5 and 3 px and reports the threshold with the best mean of
  AUC@1, AUC@3 and AUC@5; to do the same, run the benchmark once per
  threshold and keep the best.

For example, with the matches of SuperPoint and LightGlue above:

.. code-block:: bash

   python experiments/hpatches_benchmark.py data/hpatches_superpoint_lightglue.npz \
       --protocol gluefactory --threshold 0.5 --output superpoint_lightglue_gf.json


Synthetic benchmark
-------------------

What it measures
^^^^^^^^^^^^^^^^

``synthetic_benchmark.py`` solves problems from
:func:`~homlib.generate_problem_instance`. Each problem is a random plane seen
by two cameras with a random relative pose and a 70° field of view. The
generator distorts the points of each image with the division model (see
:ref:`solver_guide`), with a random parameter between -0.2 and -0.01, on the
side or sides that its :class:`~homlib.DistortionCase` sets, and returns the
points, the ground-truth homography and the distortion parameters. The benchmark adds
Gaussian noise to the points, given in pixels for a focal length of 1000 px.

The *error* of a solution is the larger of two errors, both in the
coordinates of the generator, for a focal length of 1: the homography
error, the Frobenius norm of the difference between the estimated and the
true homography, each scaled to a Frobenius norm of 1 with a positive last
element; and :math:`|k_1 - \hat k_1| + |k_2 - \hat k_2|`, the sum of the
errors of the two distortion parameters. The *best solution* is the solution
with the smallest error.

Unlike HPatches, these problems have distortion and an exact solution. The
benchmark measures each solver on the distortion cases it models (``NO_DISTORTION``
and the case it is made for, or all five for the general two-sided solvers),
with the experiments of Wadenbäck et al. [Wadenback2026]_:

**Numerical stability** (their Fig. 3)
   Without noise, how often the solver misses the true solution: the error of
   its best solution is above 1e-6. A stable solver misses almost never.

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
given shares of outliers, with a 5 px inlier threshold. The outliers replace
the points of the second image of some matches by random points; the *true
inliers* are the other matches. The experiment reports the share of the true
inliers that LOMSAC finds, and the median time of a run. Each solver
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

``--output`` and ``--compare A B`` work as in the HPatches benchmark; the saved
settings are ``--problems``, ``--noise``, ``--iterations`` and ``--outliers``.
Here the two runs don't solve the same problems, since
:func:`~homlib.generate_problem_instance` seeds itself with the time: the
differences are statistical. The differences of the rates (no solution,
missed, inliers found) come with their standard errors. The ratios of the
median errors and of the median times have none: take them as indications.
Without noise, the median errors are close to the machine precision, and
their ratio is not given (``n/a``).


Regression test
---------------

``python/tests/test_hpatches.py`` runs every point solver in LOMSAC on 8 HPatches
pairs (4 viewpoint, 4 illumination), with OpenCV SIFT matches, a 3 px inlier
threshold and 3 seeds. It checks that each solver finds, on average, at least
95% of the true inliers. The two-sided solvers of Kukelova et al. and of
Wadenbäck et al. have lower bounds, from 75% to 85%: they are weaker on the
illumination pairs, whose homography is the identity, which any equal
distortion parameters fit. The test is part of the test suite, and takes a few
seconds. To run it alone:

.. code-block:: bash

   python -m pytest tests/test_hpatches.py

Test data
^^^^^^^^^

The matches of the 8 pairs are not in the repository. They are the asset
``hpatches_sift_subset.npz`` of the GitHub release ``test-data-v1``. The test
expects the file with the SHA-256 checksum ``SUBSET_SHA256``, and looks for it,
in this order:

1. ``$HOMLIB_TEST_DATA/hpatches_sift_subset.npz``, if the environment variable
   ``HOMLIB_TEST_DATA`` names a directory that contains it. If this file has
   another checksum, the test fails;
2. a cached copy in ``~/.cache/homlib/`` (or ``$XDG_CACHE_HOME/homlib/``), if it
   has the checksum. A cached file with another checksum, such as an older
   version of the data, is downloaded again;
3. a download of ``SUBSET_URL``, the release asset, checked against the
   checksum, then cached. ``HOMLIB_TEST_DATA_URL`` replaces this URL, for
   example to test new data from a fork's release before it is published.

Without the file, the test is skipped. If ``HOMLIB_REQUIRE_TEST_DATA`` is
set, as in CI, the test fails instead, so that a missing file doesn't go
unnoticed.

Publishing the test data (maintainers)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The release must exist before CI runs the test. To write the file from the
full SIFT matches, publish it, and print its checksum:

.. code-block:: bash

   python experiments/hpatches_features.py --matches data/hpatches_sift.npz \
       --export_subset hpatches_sift_subset.npz
   gh release create test-data-v1 hpatches_sift_subset.npz \
       -R marcusvaltonen/HomLib --latest=false \
       --title "Test data" --notes "Data files for the tests (not a HomLib release)"
   shasum -a 256 hpatches_sift_subset.npz

To change the data, publish it under a new tag (``test-data-v2``), and update
``SUBSET_URL`` and ``SUBSET_SHA256`` in the test: each checkout tests the file
that its test expects, and downloads it again if the cache has another
version. HPatches is under the MIT license, which allows publishing these
matches.


Results
-------

The table gives the HPatches results of five solvers with OpenCV SIFT matches,
LOMSAC, a 3 px threshold and 3 seeds per pair, on an AWS g4dn.xlarge host
(Linux, Ubuntu 24.04): on macOS, the shares of the true inliers found differ
slightly (see ``--seeds``). The times depend on the host: compare them within
a run.

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
