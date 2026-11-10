.. _benchmarking_guide:

Benchmarking
============

HomLib's benchmarks compare its solvers and robust estimators, to choose a
solver for your data, or to measure the effect of a change in HomLib: run a
benchmark before and after the change, and compare the two runs.

.. toctree::
   :hidden:

   benchmarks/hpatches
   benchmarks/synthetic
   benchmarks/heb

- :doc:`benchmarks/hpatches`: real image pairs of planar scenes, without lens
  distortion.
- :doc:`benchmarks/synthetic`: random problems with a known solution, with or
  without distortion and noise.
- :doc:`benchmarks/heb`: real pairs of non-planar scenes from photo
  collections, where most matches are outliers.

The scripts are in ``python/benchmark/``, one directory per benchmark. They
need the ``homlib`` Python package (see :ref:`installation_users`), and
``pip install opencv-python tqdm`` (the HEB benchmark also needs ``h5py``, and
the LightGlue and LoMa matchers their own packages). Run them from the
``python/`` directory of the repository.

.. _benchmark_common_options:

Common options
--------------

``--output FILE`` and ``--compare A B``
   ``--output`` saves the settings and the results of a run. ``--compare``
   (the ``compare A B`` command for HEB) prints the differences B - A between
   two saved runs: with their standard errors for the shares and the AUCs or
   mAAs (a difference smaller than about twice its standard error could be
   chance), and as ratios for the times and median errors. It warns if the two
   runs have different settings. The times depend on the host and its load:
   compare them within a run, or between runs on the same host.

``--backend NAME`` (HPatches, HEB, and synthetic with ``--outliers``)
   The robust estimator: ``lomsac`` (default, see :ref:`backend_guide`),
   ``gcransac``, ``magsac`` or ``superansac``. The last three need HomLib built
   with them (see :ref:`backend_availability`). MAGSAC++ has no inlier
   threshold: it gets the threshold as its largest noise level,
   ``sigma_max``.

``--seeds N`` (HPatches, default 3; HEB, default 1)
   The number of runs per pair. LOMSAC gets fixed seeds, so its runs are
   reproducible on the same platform (other platforms draw other samples from
   the same seed).

Regression test
---------------

``python/tests/test_hpatches.py`` runs every point solver in LOMSAC on 8
HPatches pairs, and checks that each finds at least 95% of the true inliers
(75% to 85% for the two-sided solvers, which are weaker on the pairs whose
homography is the identity). It is part of the test suite, downloads its data
file on the first run, and is skipped without it:

.. code-block:: bash

   python -m pytest tests/test_hpatches.py
