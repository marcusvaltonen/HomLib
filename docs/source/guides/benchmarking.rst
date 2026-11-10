.. _benchmarking_guide:

Benchmarking
============

HomLib's benchmarks compare its solvers and robust estimators, to choose a
solver for your data, or to measure the effect of a change in HomLib: run a
benchmark before and after the change, and compare the two runs.

.. toctree::
   :maxdepth: 1

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
need the ``homlib`` Python package (see :ref:`installation_users`) and a few
packages of their own, and run from the ``python/`` directory of the
repository:

.. code-block:: bash

   pip install opencv-python tqdm

Common options
--------------

``--backend NAME``
   The robust estimator: ``lomsac`` (default, see :ref:`backend_guide`),
   ``gcransac``, ``magsac`` or ``superansac``. The last three need HomLib built
   with them (see :ref:`backend_availability`). MAGSAC++ has no inlier
   threshold: it gets the threshold as its largest noise level,
   ``sigma_max``.

``--output FILE`` and ``--compare A B``
   ``--output`` saves the settings and the results of a run. ``--compare``
   prints the differences B - A between two saved runs, with their standard
   errors: a difference smaller than about twice its standard error could be
   chance. It warns if the two runs have different settings.

``--seeds N``
   The number of runs per problem or pair. LOMSAC gets fixed seeds, so its runs
   are reproducible on the same platform (other platforms draw other samples
   from the same seed).

The times depend on the host and its load: compare them within a run.

Regression test
---------------

``python/tests/test_hpatches.py`` runs every point solver in LOMSAC on 8
HPatches pairs, and checks that each finds at least 95% of the true inliers
(75% to 85% for the two-sided solvers, which are weaker on the pairs whose
homography is the identity). It is part of the test suite:

.. code-block:: bash

   python -m pytest tests/test_hpatches.py

Test data
^^^^^^^^^

The matches of the 8 pairs are the asset ``hpatches_sift_subset.npz`` of the
GitHub release ``test-data-v1``, checked against the checksum
``SUBSET_SHA256`` of the test. The test uses
``$HOMLIB_TEST_DATA/hpatches_sift_subset.npz`` if this variable is set, or a
copy cached in ``~/.cache/homlib/``, or downloads it (from
``HOMLIB_TEST_DATA_URL`` if set). Without the file, the test is skipped, or
fails if ``HOMLIB_REQUIRE_TEST_DATA`` is set, as in CI.

To change the data, maintainers write the new file, publish it under a new
release tag, and update ``SUBSET_URL`` and ``SUBSET_SHA256`` in the test:

.. code-block:: bash

   python benchmark/hpatches/hpatches_features.py --matches data/hpatches_sift.npz \
       --export_subset hpatches_sift_subset.npz
   gh release create test-data-v2 hpatches_sift_subset.npz -R marcusvaltonen/HomLib \
       --latest=false --title "Test data" --notes "Data files for the tests (not a HomLib release)"
   shasum -a 256 hpatches_sift_subset.npz
