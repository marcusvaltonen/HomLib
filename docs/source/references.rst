.. _references:

References
==========

If you use HomLib in your work, please cite the relevant paper(s) below.

Point-based solvers
-------------------

- :class:`~homlib.CaiPAMI2025` :cite:`Cai2025`
- :class:`~homlib.DLT` :cite:`HartleyZisserman2003`
- :class:`~homlib.GuoPRL2019` :cite:`Guo2019`
- :class:`~homlib.FitzgibbonCVPR2001OneSided`,
  :class:`~homlib.FitzgibbonCVPR2001TwoSidedEqual` :cite:`Fitzgibbon2001`
- :class:`~homlib.KukelovaCVPR2015TwoSided`,
  :class:`~homlib.KukelovaCVPR2015TwoSided6Pt`,
  :class:`~homlib.KukelovaCVPR2015TwoSidedEqual`,
  :class:`~homlib.KukelovaCVPR2015TwoSidedEqual6Pt` :cite:`Kukelova2015`
- :cite:`ValtonenOrnhag2020`
- :cite:`ValtonenOrnhag2021`
- :class:`~homlib.NakanoICPR2025OneSided`,
  :class:`~homlib.NakanoICPR2025OneSidedRight` :cite:`Nakano2024`.
  Note: The class name uses ``ICPR2025`` (the conference year), while the
  publication is dated 2024.
- :cite:`Barath2023heb`
- :class:`~homlib.Wadenback3DV2026OneSided`,
  :class:`~homlib.Wadenback3DV2026TwoSidedEqual`,
  :class:`~homlib.Wadenback3DV2026TwoSided` :cite:`Wadenback2026`

Affine-covariant and orientation-covariant solvers
---------------------------------------------------

- :class:`~homlib.BarathVISAPP2016Affine` :cite:`Barath2016`
- :class:`~homlib.ValtonenOrnhagICPR2026AffineOneSided`,
  :class:`~homlib.ValtonenOrnhagICPR2026OrientationOneSided`
  :cite:`ValtonenOrnhagAdalbjornsson2026`

Robust estimation
-----------------

- :class:`~homlib.LOMSAC` :cite:`Lebeda2012`
- :class:`~homlib.GCRANSAC` :cite:`Barath2018`
- :cite:`Barath2019magsac`
- :class:`~homlib.MAGSAC` :cite:`Barath2020magsacpp`
- :class:`~homlib.SupeRANSAC` :cite:`Barath2025superansac`

Automatic solver generation
----------------------------

- :cite:`Larsson2017`

Miscellaneous
-------------

- :cite:`LopezAntequera2019` — The distortion of real cameras, for the range
  check on the distortion (:meth:`~homlib.Solver.set_k_range`).

.. note::

   OpticalLimits lens reviews (distortion without in-camera correction):
   `Nikkor Z 70-200mm f/2.8 VR S II
   <https://opticallimits.com/nikon/nikon-z/nikkor-z-70-200mm-f-2-8-vr-s-ii-review/>`_ and
   `Sigma 20-200mm f/3.5-6.3 DG Contemporary
   <https://opticallimits.com/sony/sony-fe-full-format/review-of-the-sigma-20-200mm-f-3-5-6-3-dg-contemporary-in-sony-e-mount/>`_.
   Typical distortion of long lenses, for the range check on the distortion.

.. rubric:: Full bibliography

.. bibliography::
   :all:
