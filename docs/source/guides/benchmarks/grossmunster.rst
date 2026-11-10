Grossmünster Church benchmark
=============================

Wadenbäck et al. [Wadenback2026]_ (Section 5, Table 3) evaluate the two-sided
solvers on six pairs of fisheye images of the Grossmünster Church in Zürich,
from `RadialSfM <https://github.com/vlarsson/radialsfm>`__ (Larsson et al.,
ECCV 2020). Unlike HPatches, both images of each pair have a strong
distortion, which the division model only approximates. The scripts in
``python/benchmark/grossmunster/`` rebuild the paper's ground truth, run its
seven two-sided solvers, and draw its Figures 7 and 8.

The ground truth of a pair is the homography of the wall with the two towers,
:math:`H = K (R + t n^\top / d) K^{-1}`, from the relative pose
:math:`(R, t)` of the two cameras and the plane :math:`n^\top x = d` of the
wall in the first camera. ``grossmunster_ground_truth.py`` gets them from a 3D
reconstruction of the images DSC_2510 to DSC_2522 (pycolmap), with the camera
of the paper's authors (a fisheye camera, fixed), and a plane fitted to the
wall's 3D points.

The benchmark runs each solver in :class:`~homlib.LOMSAC` with the paper's
settings, as its authors gave them: OpenCV SIFT matches (8000 keypoints, ratio
test at 0.8), a 5 px inlier threshold, at least 500 iterations, and a final
least-squares refinement. It reports, as means over the pairs and 10 runs per
pair:

**Error**
   The paper's error :math:`\|H - \hat H\|_F`, with both homographies scaled
   to a unit norm and a positive determinant, in pixel coordinates.
   ``Error (cal.)`` measures it in calibrated coordinates
   (:math:`K^{-1} H K`), where the translation doesn't dominate the norm.

**Inliers**, **Time**
   The number of inliers found, and the time of a run of the estimator.

The pairs are the five that the authors remember, DSC_2512/2513, 2513/2514,
2514/2515, 2517/2518 and 2518/2519, and 2519/2520, the likely sixth.
``--pairs`` selects others.

Running it
----------

Download the images of RadialSfM's "Grossmunster Church" from the Google Drive
link of `RadialSfM's README <https://github.com/vlarsson/radialsfm>`__ (Google
Drive refuses scripted downloads), and unzip them in ``data/``. Then, with
``pip install pycolmap matplotlib`` for the ground truth and the figures:

.. code-block:: bash

   python benchmark/grossmunster/grossmunster_ground_truth.py --data_dir data/fisheye_grossmunster \
       --sfm_dir data/grossmunster_sfm --output data/grossmunster_sift.npz
   python benchmark/grossmunster/grossmunster_benchmark.py data/grossmunster_sift.npz --per_pair
   python benchmark/grossmunster/grossmunster_figures.py data/grossmunster_sift.npz \
       --data_dir data/fisheye_grossmunster --output_dir figures

The first command takes a few minutes. It prints the reprojection errors of
the reconstruction (0.9 px on average), and checks the ground truth on the
reconstruction's points of the wall. On the six pairs, their median distance
to the ground truth is 0.7 to 1.7 px, against 0.65 to 1.3 px for the best
homography fitted to them: the wall isn't flat, so no homography fits it
exactly. ``--intrinsics refined`` calibrates the camera in the reconstruction
instead of using the authors' one. This halves the reprojection errors (0.5 px
on average) and brings the ground truth within 0.1 px of the best homography,
but changes the benchmark's errors by less than 0.007.

.. _grossmunster_results:

Results
-------

The seven solvers of Table 3 with the default settings, on a MacBook Pro (M5
Pro), with the commands above:

.. list-table::
   :header-rows: 1

   * - Solver
     - Error
     - Error (cal.)
     - Time (ms)
     - Paper: error
     - Paper: time (ms)
   * - :class:`~homlib.KukelovaCVPR2015TwoSided`
     - 0.014
     - 0.023
     - 26
     - 0.039
     - 2657
   * - :class:`~homlib.KukelovaCVPR2015TwoSided6Pt`
     - 0.014
     - 0.023
     - 37
     - 0.040
     - 330
   * - :class:`~homlib.Wadenback3DV2026TwoSided`
     - 0.014
     - 0.022
     - 15
     - 0.038
     - 150
   * - :class:`~homlib.FitzgibbonCVPR2001TwoSidedEqual`
     - 0.011
     - 0.009
     - 79
     - 0.040
     - 233
   * - :class:`~homlib.KukelovaCVPR2015TwoSidedEqual`
     - 0.012
     - 0.014
     - 42
     - 0.040
     - 2053
   * - :class:`~homlib.KukelovaCVPR2015TwoSidedEqual6Pt`
     - 0.014
     - 0.015
     - 54
     - 0.039
     - 321
   * - :class:`~homlib.Wadenback3DV2026TwoSidedEqual`
     - 0.010
     - 0.009
     - 11
     - 0.039
     - 84

As in the paper, the solvers' errors are close to each other, and Wadenbäck
et al.'s solvers are the fastest of each kind. The errors are about a third
of the paper's: most of them come from one pair, 2519/2520 (``--per_pair``),
which may not be the paper's sixth pair, and the paper's ground truth may
differ from this one. The times are those of a different computer.
