.. _references:

References
==========

If you use HomLib in your work, please cite the relevant paper(s) below.

Point-based solvers
-------------------

.. [Fitzgibbon2001] A. Fitzgibbon, "Simultaneous linear estimation of multiple view geometry
   and lens distortion," in *Proceedings of the IEEE Computer Society Conference on Computer
   Vision and Pattern Recognition (CVPR)*, 2001.

   See :func:`~homlib.estimate_fitzgibbon_cvpr_2001_one_sided`,
   :func:`~homlib.estimate_fitzgibbon_cvpr_2001_two_sided_equal`.

.. [Kukelova2015] Z. Kukelova, J. Heller, M. Bujnak, and T. Pajdla, "Radial Distortion
   Homography," in *Proceedings of the IEEE Conference on Computer Vision and Pattern
   Recognition (CVPR)*, 2015.

   See :func:`~homlib.estimate_kukelova_cvpr_2015_two_sided`,
   :func:`~homlib.estimate_kukelova_cvpr_2015_two_sided_6pt`.

.. [ValtonenOrnhag2020] M. Valtonen Örnhag, P. Persson, M. Wadenbäck, K. Åström, and
   A. Heyden, "Minimal Solvers for Indoor UAV Positioning," in *Proceedings of the
   International Conference on Pattern Recognition (ICPR)*, 2021, pp. 1136--1143.

.. [ValtonenOrnhag2021] M. Valtonen Örnhag, P. Persson, M. Wadenbäck, K. Åström, and
   A. Heyden, "Efficient Real-Time Radial Distortion Correction for UAVs," in *Proceedings
   of the IEEE/CVF Winter Conference on Applications of Computer Vision (WACV)*, 2021,
   pp. 1751--1760.

.. [Nakano2024] G. Nakano, "Inverse DLT Method for One-Sided Radial Distortion Homography,"
   in *Proceedings of the International Conference on Pattern Recognition (ICPR)*, 2024,
   pp. 448--462.

   See :func:`~homlib.estimate_nakano_icpr_2025_one_sided`.
   Note: The function name uses ``icpr_2025`` (the conference year), while the
   publication is dated 2024.

.. [Wadenback2026] M. Wadenbäck, M. Valtonen Örnhag, and J. Edstedt, "Radially Distorted
   Homographies, Revisited," in *Proceedings of the International Conference on 3D Vision
   (3DV)*, 2026.

   See :func:`~homlib.estimate_wadenback_3dv_2026_one_sided`,
   :func:`~homlib.estimate_wadenback_3dv_2026_two_sided_equal`,
   :func:`~homlib.estimate_wadenback_3dv_2026_two_sided`.

Affine-covariant and orientation-covariant solvers
---------------------------------------------------

.. [Barath2016] D. Barath and L. Hajder, "Novel Ways to Estimate Homography from Local Affine
   Transformations," in *Proceedings of the International Conference on Computer Vision Theory
   and Applications (VISAPP)*, 2016.

   See :func:`~homlib.lomsac_barath_visapp_2016_affine`.

.. [ValtonenOrnhagAdalbjornsson2026] M. Valtonen Örnhag and S. Adalbjörnsson, "Radial
   Distortion Homography Estimation from Affine-Covariant or Orientation-Covariant Features,"
   in *Proceedings of the International Conference on Pattern Recognition (ICPR)*, 2026.

   See :func:`~homlib.lomsac_valtonenornhag_icpr_2026_one_sided_affine`,
   :func:`~homlib.lomsac_valtonenornhag_icpr_2026_one_sided_ori`.

Robust estimation
-----------------

.. [Lebeda2012] K. Lebeda, J. Matas, and O. Chum, "Fixing the Locally Optimized RANSAC," in
   *Proceedings of the British Machine Vision Conference (BMVC)*, 2012.

Automatic solver generation
----------------------------

.. [Larsson2017] V. Larsson, K. Åström, and M. Oskarsson, "Efficient Solvers for Minimal
   Problems by Syzygy-based Reduction," in *Proceedings of the IEEE Conference on Computer
   Vision and Pattern Recognition (CVPR)*, 2017.
