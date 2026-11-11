HPatches benchmark
==================

`HPatches <https://github.com/hpatches/hpatches-dataset>`__ (Balntas et al.,
CVPR 2017) has 116 sequences of 6 images of planar scenes: in the 59
*viewpoint* sequences the camera moves, in the 57 *illumination* sequences the
lighting changes. Each sequence gives the ground-truth homography from its
first image to the 5 others: 580 pairs. The images have no lens distortion, so
the benchmark measures how each solver does on real matches, and what the
distortion solvers lose by estimating parameters that are zero.

It runs a robust estimator (LOMSAC by default) with every point solver on every
pair, 3 times per pair, and reports, for the viewpoint and the illumination
pairs:

**Inliers found**
   The share of the true inliers (the matches within the inlier threshold, 3 px
   by default, of the ground-truth homography) that the estimator returns.
   It ignores the false inliers: the corner error accounts for them. The pairs
   without any true inlier are left out (except with ``--protocol
   gluefactory``).

**AUC@1, AUC@3, AUC@5, AUC@10**
   The corner error of a pair maps the four corners of the first image with
   the estimated model (distortion included) and with the ground truth, and
   averages the distances. The AUC@t is the area under its cumulative curve up
   to t pixels, in percent: 100 for no error, 0 for an error above t or a
   failed estimate. Wadenbäck et al. [Wadenback2026]_ report it as "mAA (AUC)".

**Time**
   The median time of a run of the estimator.

Running it
----------

The benchmark reads a matches file: the matches of the 580 pairs, a score per
match, and the ground truth, so that every run compares the solvers on the
same matches. The release
`test-data-v1 <https://github.com/marcusvaltonen/HomLib/releases/tag/test-data-v1>`__
of HomLib has three of them, each with its SHA-256:

- ``hpatches_sift.npz``: OpenCV SIFT, with the defaults of
  ``hpatches_features.py`` (below);
- ``hpatches_superpoint_lightglue.npz``: SuperPoint and LightGlue with the
  setup of Wadenbäck et al. (below);
- ``hpatches_loma.npz``: LoMa, 2048 keypoints.

The score of a match is higher for a better match: 1 minus the distance ratio
of the ratio test for SIFT, and the match confidence for LightGlue and LoMa.
The benchmark doesn't use them yet: they are for the ``scores`` argument of
the robust estimators' ``estimate``, which the samplers guided by the quality
of the matches use.

.. code-block:: bash

   mkdir -p data && cd data
   for f in hpatches_sift hpatches_superpoint_lightglue hpatches_loma; do
       curl -LO https://github.com/marcusvaltonen/HomLib/releases/download/test-data-v1/$f.npz
       curl -LO https://github.com/marcusvaltonen/HomLib/releases/download/test-data-v1/$f.npz.sha256
       shasum -a 256 -c $f.npz.sha256
   done
   cd ..

To compute matches of your own, ``hpatches_features.py`` downloads the
HPatches images (about 2 GB) if needed, and writes a matches file:

.. code-block:: bash

   python benchmark/hpatches/hpatches_features.py --data_dir data --output data/hpatches_sift.npz

The default is OpenCV SIFT (4000 keypoints, ratio test at 0.8). ``--matcher
lightglue`` with ``--features superpoint``, ``disk``, ``aliked`` or ``sift``
uses LightGlue, and ``--features loma`` uses LoMa (see
:ref:`their installation <benchmark_installation>`).
The setup of Wadenbäck et al. is SuperPoint and LightGlue with 1024 keypoints
on images resized to 480 px:

.. code-block:: bash

   python benchmark/hpatches/hpatches_features.py --data_dir data --features superpoint \
       --matcher lightglue --num_features 1024 --resize 480 \
       --output data/hpatches_superpoint_lightglue.npz

Then run the solvers, before and after a change:

.. code-block:: bash

   python benchmark/hpatches/hpatches_benchmark.py data/hpatches_sift.npz --output before.json
   # ... change HomLib, and reinstall it ...
   python benchmark/hpatches/hpatches_benchmark.py data/hpatches_sift.npz --output after.json
   python benchmark/hpatches/hpatches_benchmark.py --compare before.json after.json

Besides the :ref:`common options <benchmark_common_options>`, ``--solvers`` selects
the solvers, ``--threshold`` sets the inlier threshold in pixels, and
``--iterations N`` runs exactly N iterations, so that every solver gets the same
budget (as Fig. 6 of Wadenbäck et al.). ``--protocol gluefactory`` follows
`glue-factory <https://github.com/cvg/glue-factory>`__'s HPatches evaluation,
so that the AUCs compare with its results: images resized to a shorter side of
480 px, its ground-truth convention, 540 pairs, and every pair counted.
``--help`` lists all the options, and the scripts' docstrings give the details.

.. _hpatches_results:

Results
-------

Five of the solvers, with OpenCV SIFT matches, LOMSAC, a 3 px threshold and 3
seeds per pair, on an AWS g4dn.xlarge host (Ubuntu 24.04), with the default
settings of both scripts:

.. code-block:: bash

   python benchmark/hpatches/hpatches_features.py --data_dir data --output data/hpatches_sift.npz
   python benchmark/hpatches/hpatches_benchmark.py data/hpatches_sift.npz

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
accurate, and the more distortion parameters a solver estimates, the lower its
AUC. The two-sided solver loses most on the illumination pairs, whose
homography is the identity, which any equal distortion parameters fit.
