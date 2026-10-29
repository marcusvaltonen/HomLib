.. _api:

homlib API
=================

.. currentmodule:: homlib

Data types
----------

.. autosummary::
   :toctree: generated/

    PoseData
    LORansacOptions
    RansacStatistics
    DistortionCase
    ProblemConfig
    ProblemInstance

Point-based solvers (minimal)
-----------------------------

Direct solver calls that return all algebraic solutions from a minimal sample.
Distortion suffixes: ``one_sided`` — single distortion coefficient on the left image,
``two_sided_equal`` — equal distortion on both images, ``two_sided`` — independent distortion on each image.

.. autosummary::
   :toctree: generated/

    estimate_fitzgibbon_cvpr_2001_one_sided
    estimate_fitzgibbon_cvpr_2001_two_sided_equal
    estimate_kukelova_cvpr_2015_two_sided
    estimate_kukelova_cvpr_2015_two_sided_6pt
    estimate_nakano_icpr_2025_one_sided
    estimate_wadenback_3dv_2026_one_sided
    estimate_wadenback_3dv_2026_two_sided_equal
    estimate_wadenback_3dv_2026_two_sided

Point-based solvers (LOMSAC)
----------------------------

Robust estimation using Locally Optimized MSAC (LOMSAC).
These functions take a full set of correspondences and return a single best estimate
together with inlier statistics.

.. autosummary::
   :toctree: generated/

    lomsac_fitzgibbon_cvpr_2001_one_sided
    lomsac_fitzgibbon_cvpr_2001_two_sided_equal
    lomsac_kukelova_cvpr_2015_two_sided
    lomsac_kukelova_cvpr_2015_two_sided_equal
    lomsac_kukelova_cvpr_2015_two_sided_6pt
    lomsac_kukelova_cvpr_2015_two_sided_equal_6pt
    lomsac_nakano_icpr_2025_one_sided
    lomsac_wadenback_3dv_2026_one_sided
    lomsac_wadenback_3dv_2026_two_sided_equal
    lomsac_wadenback_3dv_2026_two_sided

Affine and orientation solvers (LOMSAC)
---------------------------------------

Solvers that leverage affine-covariant (A) or orientation-covariant (O) features,
requiring fewer correspondences than point-based solvers. Only LOMSAC wrappers are
provided since these are non-minimal solvers.

.. autosummary::
   :toctree: generated/

    lomsac_barath_visapp_2016_affine
    lomsac_valtonenornhag_icpr_2026_one_sided_affine
    lomsac_valtonenornhag_icpr_2026_one_sided_ori

Synthetic data generation
-------------------------

Utilities for generating synthetic problem instances with known ground truth,
useful for benchmarking and testing.

.. autosummary::
   :toctree: generated/

    generate_problem_instance
