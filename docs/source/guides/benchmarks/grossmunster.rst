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
wall in the first camera. They come from a 3D reconstruction of the whole
scene (its 373 images, with pycolmap's global mapper), which calibrates the
fisheye camera and gives the poses; the plane is fitted to the 3D points that the
images of the wall, DSC_2510 to DSC_2522, see.

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

The release
`test-data-v1 <https://github.com/marcusvaltonen/HomLib/releases/tag/test-data-v1>`__
of HomLib has the benchmark's data, each file with its SHA-256:

- ``grossmunster_sift.npz`` (0.6 MB): the SIFT matches and the ground truth
  of the pairs, all that the benchmark needs;
- ``grossmunster_sfm.zip`` (27 MB): the reconstruction of the whole scene, in
  COLMAP's format, to recompute the ground truth or for other uses.

.. code-block:: bash

   cd data
   for f in grossmunster_sift.npz grossmunster_sfm.zip; do
       curl -LO https://github.com/marcusvaltonen/HomLib/releases/download/test-data-v1/$f
       curl -LO https://github.com/marcusvaltonen/HomLib/releases/download/test-data-v1/$f.sha256
       shasum -a 256 -c $f.sha256
   done
   unzip grossmunster_sfm.zip
   cd ..
   python benchmark/grossmunster/grossmunster_benchmark.py data/grossmunster_sift.npz --per_pair

The ground truth and the figures need the images: download RadialSfM's
"Grossmunster Church" from the Google Drive link of
`RadialSfM's README <https://github.com/vlarsson/radialsfm>`__ (Google Drive
refuses scripted downloads), unzip it in ``data/``, and
``pip install pycolmap matplotlib``:

.. code-block:: bash

   python benchmark/grossmunster/grossmunster_ground_truth.py --data_dir data/fisheye_grossmunster \
       --sfm_model data/grossmunster_sfm --output data/grossmunster_sift.npz
   python benchmark/grossmunster/grossmunster_figures.py data/grossmunster_sift.npz \
       --data_dir data/fisheye_grossmunster --output_dir figures

The ground truth script prints the reprojection errors of the images of the
wall (0.5 px on average), and checks the ground truth on the reconstruction's
points of the wall. On the six pairs, their median distance to the ground
truth is 0.6 to 1.0 px, within 0.1 px of the best homography fitted to them:
the wall isn't flat, so no homography fits it exactly, but the ground truth is
about as good as one can be. ``grossmunster_scene_sfm.py`` recomputes the
reconstruction of the whole scene (about 20 minutes). Without
``--sfm_model``, the ground truth script reconstructs the images of the wall
only (``--sfm_dir``), with the authors' camera (``--intrinsics paper``) or a
refined one: both change the benchmark's errors by less than 0.006.

.. _grossmunster_results:

Results
-------

The seven solvers of Table 3 with the default settings, on a MacBook Pro (M5
Pro), with the downloaded data and the command above:

.. list-table::
   :header-rows: 1

   * - Solver
     - Error
     - Error (cal.)
     - Time (ms)
     - Paper: error
     - Paper: time (ms)
   * - :class:`~homlib.KukelovaCVPR2015TwoSided`
     - 0.017
     - 0.024
     - 26
     - 0.039
     - 2657
   * - :class:`~homlib.KukelovaCVPR2015TwoSided6Pt`
     - 0.015
     - 0.022
     - 37
     - 0.040
     - 330
   * - :class:`~homlib.Wadenback3DV2026TwoSided`
     - 0.016
     - 0.023
     - 15
     - 0.038
     - 150
   * - :class:`~homlib.FitzgibbonCVPR2001TwoSidedEqual`
     - 0.011
     - 0.010
     - 78
     - 0.040
     - 233
   * - :class:`~homlib.KukelovaCVPR2015TwoSidedEqual`
     - 0.013
     - 0.012
     - 42
     - 0.040
     - 2053
   * - :class:`~homlib.KukelovaCVPR2015TwoSidedEqual6Pt`
     - 0.009
     - 0.014
     - 56
     - 0.039
     - 321
   * - :class:`~homlib.Wadenback3DV2026TwoSidedEqual`
     - 0.010
     - 0.010
     - 12
     - 0.039
     - 84

As in the paper, the solvers' errors are close to each other, and Wadenbäck
et al.'s solvers are the fastest of each kind. The errors are a quarter to a half
of the paper's: most of them come from one pair, 2519/2520 (``--per_pair``),
which may not be the paper's sixth pair, and the paper's ground truth may
differ from this one. The times are those of a different computer.
