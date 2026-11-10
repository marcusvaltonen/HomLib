HEB benchmark
=============

HEB, the Homography Estimation Benchmark of Barath, Mishkin et al.
[Barath2023]_, has image pairs of 10 scenes of photo collections, with their
SIFT matches, the matches that are inliers of the ground-truth homography, and
the cameras. The scenes aren't planar: each homography maps one plane of the
scene, so most matches are outliers.

``heb_benchmark.py`` follows the paper's protocol (its Section 4). For each
homography, it keeps the matches whose SNN ratio (Lowe's ratio: the descriptor
distance to the nearest neighbour over the distance to the second nearest) is
below a threshold, runs a robust estimator, and measures mAAs (mean average
accuracies: for each of 10 error thresholds, the share of the homographies with
an error below it, averaged over the thresholds), per scene, then averaged over
the scenes:

**mAA (10°, 5 m)**, the main metric
   The average of the mAAs of the rotation error (1 to 10°) and of the
   absolute translation error (0.1 to 5 m) of the pose decomposed from the
   homography.

**mAA of the transfer error**
   The mean distance, in image 2, between the ground-truth inliers of image 1
   mapped by the estimated model and their matches (1 to 20 px; the
   benchmark's code calls it the reprojection error).

**Time**
   The mean time of a run of the estimator.

Running it
----------

The dataset (33 GB, of which 25 GB for the test split) is linked from the
`benchmark's README <https://github.com/danini/homography-benchmark>`__: put it
in a directory with ``test/`` and ``training_and_validation/`` (or ``train/``). The protocol
has two steps: ``tune`` searches, for each solver, the best inlier threshold
and SNN threshold on the 2 training scenes, with 1000 iterations; ``test``
runs each solver with these settings on the 8 test scenes, for exactly 10,
100, 1000 and 10,000 iterations (``--iterations``):

.. code-block:: bash

   pip install h5py
   python benchmark/heb/heb_benchmark.py tune /path/to/heb --output settings.json
   python benchmark/heb/heb_benchmark.py test /path/to/heb --settings settings.json --output heb.json
   python benchmark/heb/heb_benchmark.py compare before.json after.json

The full dataset takes days on 4 cores, so ``--max_pairs`` evaluates a fixed
random subset of the homographies of each scene (250 for ``tune``, 1000 for
``test``; 0 for all, as the paper). Without ``--settings``, ``test`` uses a
4 px threshold and an SNN threshold of 0.8 for every solver; ``--snn 1`` tests
without prefiltering.

Unlike the paper, only the two thresholds are tuned, and LOMSAC samples
uniformly, so the results compare with the paper's curves without PROSAC. The
benchmark's code has no license, so the script reimplements its metrics; on the
same estimates (without distortion), both give the same errors. The script's
docstring lists the other details.

.. _heb_results:

Results
-------

A smoke run on four solvers, on part of the data (``tune`` on 50 homographies
per training scene, ``test`` on 100 per test scene), with LOMSAC, on an AWS
g4dn.xlarge host (Linux, 4 processes in parallel). Each row gives the solver's
tuned thresholds, its mAA (10°, 5 m) in percent for each number of iterations,
and the mean time of a run at 1000 iterations:

.. list-table::
   :header-rows: 1

   * - Solver (threshold, SNN)
     - 10 it.
     - 100 it.
     - 1000 it.
     - 10,000 it.
     - Time (1000 it.)
   * - :class:`~homlib.DLT` (6 px, 0.5)
     - 27.1
     - 30.6
     - 33.4
     - 34.2
     - 16.1 ms
   * - :class:`~homlib.GuoPRL2019` (6 px, 0.5)
     - 27.1
     - 30.6
     - 33.4
     - 34.2
     - 5.5 ms
   * - :class:`~homlib.Wadenback3DV2026TwoSidedEqual` (4 px, 0.5)
     - 21.3
     - 27.3
     - 29.8
     - 31.3
     - 11.0 ms
   * - :class:`~homlib.Wadenback3DV2026TwoSided` (3 px, 0.5)
     - 17.4
     - 25.0
     - 28.7
     - 29.6
     - 19.4 ms

The accuracy grows with the number of iterations. On these scenes, which have
little distortion, the solvers without distortion score best. DLT and Guo score
the same, because a homography through 4 points is unique: from the same
samples, both give the same homography; Guo is faster.
