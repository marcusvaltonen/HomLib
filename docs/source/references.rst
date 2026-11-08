.. _references:

References
==========

If you use HomLib in your work, please cite the relevant paper(s) below.

Point-based solvers
-------------------

.. [Cai2025] S. Cai, Z. Wu, L. Guo, J. Wang, S. Zhang, J. Yan, and S. Shen, "Fast and
   Interpretable 2D Homography Decomposition: Similarity-Kernel-Similarity and
   Affine-Core-Affine Transformations," *IEEE Transactions on Pattern Analysis and Machine
   Intelligence (PAMI)*, 2025.

   See :class:`~homlib.CaiPAMI2025`.

.. [HartleyZisserman2003] R. Hartley and A. Zisserman, *Multiple View Geometry in Computer Vision*,
   2nd ed., Cambridge University Press, 2003.

   See :class:`~homlib.DLT`.

.. [Guo2019] Y. Guo, "A Minimal Solution for Plane Motion and Structure from Two Perspective
   Views," *Pattern Recognition Letters*, vol. 120, pp. 96--103, 2019.

   See :class:`~homlib.GuoPRL2019`.

.. [Fitzgibbon2001] A. Fitzgibbon, "Simultaneous linear estimation of multiple view geometry
   and lens distortion," in *Proceedings of the IEEE Computer Society Conference on Computer
   Vision and Pattern Recognition (CVPR)*, 2001.

   See :class:`~homlib.FitzgibbonCVPR2001OneSided`,
   :class:`~homlib.FitzgibbonCVPR2001TwoSidedEqual`.

.. [Kukelova2015] Z. Kukelova, J. Heller, M. Bujnak, and T. Pajdla, "Radial Distortion
   Homography," in *Proceedings of the IEEE Conference on Computer Vision and Pattern
   Recognition (CVPR)*, 2015.

   See :class:`~homlib.KukelovaCVPR2015TwoSided`,
   :class:`~homlib.KukelovaCVPR2015TwoSided6Pt`,
   :class:`~homlib.KukelovaCVPR2015TwoSidedEqual`,
   :class:`~homlib.KukelovaCVPR2015TwoSidedEqual6Pt`.

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

   See :class:`~homlib.NakanoICPR2025OneSided`,
   :class:`~homlib.NakanoICPR2025OneSidedRight`.
   Note: The class name uses ``ICPR2025`` (the conference year), while the
   publication is dated 2024.

.. [Wadenback2026] M. Wadenbäck, M. Valtonen Örnhag, and J. Edstedt, "Radially Distorted
   Homographies, Revisited," in *Proceedings of the International Conference on 3D Vision
   (3DV)*, 2026.

   See :class:`~homlib.Wadenback3DV2026OneSided`,
   :class:`~homlib.Wadenback3DV2026TwoSidedEqual`,
   :class:`~homlib.Wadenback3DV2026TwoSided`.

Affine-covariant and orientation-covariant solvers
---------------------------------------------------

.. [Barath2016] D. Barath and L. Hajder, "Novel Ways to Estimate Homography from Local Affine
   Transformations," in *Proceedings of the International Conference on Computer Vision Theory
   and Applications (VISAPP)*, 2016.

   See :class:`~homlib.BarathVISAPP2016Affine`.

.. [ValtonenOrnhagAdalbjornsson2026] M. Valtonen Örnhag and S. Adalbjörnsson, "Radial
   Distortion Homography Estimation from Affine-Covariant or Orientation-Covariant Features,"
   in *Proceedings of the International Conference on Pattern Recognition (ICPR)*, 2026.

   See :class:`~homlib.ValtonenOrnhagICPR2026AffineOneSided`,
   :class:`~homlib.ValtonenOrnhagICPR2026OrientationOneSided`.

Robust estimation
-----------------

.. [Lebeda2012] K. Lebeda, J. Matas, and O. Chum, "Fixing the Locally Optimized RANSAC," in
   *Proceedings of the British Machine Vision Conference (BMVC)*, 2012.

   See :class:`~homlib.LOMSAC`.

.. [Barath2018] D. Barath and J. Matas, "Graph-Cut RANSAC," in *Proceedings of the
   IEEE Conference on Computer Vision and Pattern Recognition (CVPR)*, 2018.

   See :class:`~homlib.GCRANSAC`.

.. [Barath2019] D. Barath, J. Matas, and J. Noskova, "MAGSAC: Marginalizing Sample
   Consensus," in *Proceedings of the IEEE Conference on Computer Vision and Pattern
   Recognition (CVPR)*, 2019.

.. [Barath2020] D. Barath, J. Noskova, M. Ivashechkin, and J. Matas, "MAGSAC++, a Fast,
   Reliable and Accurate Robust Estimator," in *Proceedings of the IEEE Conference on
   Computer Vision and Pattern Recognition (CVPR)*, 2020.

   See :class:`~homlib.MAGSAC`.

.. [Barath2025] D. Barath, "SupeRANSAC: One RANSAC to Rule Them All,"
   *arXiv preprint arXiv:2506.04803*, 2025.

   See :class:`~homlib.SupeRANSAC`.
   `arXiv <https://arxiv.org/abs/2506.04803>`__

Automatic solver generation
----------------------------

.. [Larsson2017] V. Larsson, K. Åström, and M. Oskarsson, "Efficient Solvers for Minimal
   Problems by Syzygy-based Reduction," in *Proceedings of the IEEE Conference on Computer
   Vision and Pattern Recognition (CVPR)*, 2017.

Miscellaneous
-------------

.. [LopezAntequera2019] M. López-Antequera, R. Marí, P. Gargallo, Y. Kuang, J. Gonzalez-Jimenez,
   and G. Haro, "Deep Single Image Camera Calibration with Radial Distortion," in *Proceedings of
   the IEEE/CVF Conference on Computer Vision and Pattern Recognition (CVPR)*, 2019,
   pp. 11809--11817. https://research.mapillary.com/publication/cvpr19d

   The distortion of real cameras, for the range check on the distortion
   (:meth:`~homlib.Solver.set_k_range`).

.. [OpticalLimits] OpticalLimits, lens reviews (distortion without in-camera correction):
   `Nikkor Z 70-200mm f/2.8 VR S II
   <https://opticallimits.com/nikon/nikon-z/nikkor-z-70-200mm-f-2-8-vr-s-ii-review/>`_ and
   `Sigma 20-200mm f/3.5-6.3 DG Contemporary
   <https://opticallimits.com/sony/sony-fe-full-format/review-of-the-sigma-20-200mm-f-3-5-6-3-dg-contemporary-in-sony-e-mount/>`_.

   Typical distortion of long lenses, for the range check on the distortion.
