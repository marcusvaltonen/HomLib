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

The scripts are in ``python/benchmark/``, one directory per benchmark. Run
them from the ``python/`` directory of the repository.

.. _benchmark_installation:

Installation
------------

The scripts need the ``homlib`` Python package (see :ref:`installation_users`)
and a few packages (``h5py`` is for the HEB benchmark):

.. code-block:: bash

   pip install opencv-python tqdm h5py

The benchmarks that match images themselves (HPatches) use OpenCV SIFT by
default, or learned matchers: `LightGlue <https://github.com/cvg/LightGlue>`__ (with
SuperPoint, DISK, ALIKED or SIFT features) and
`LoMa <https://github.com/davnords/LoMa>`__. Both run much faster on a GPU, and
install from a clone of their repository, as their READMEs explain
(`LightGlue's <https://github.com/cvg/LightGlue#installation-and-demo-->`__,
`LoMa's <https://github.com/davnords/LoMa#setupinstall>`__); the PyPI package
named ``loma`` is another project:

.. code-block:: bash

   git clone https://github.com/cvg/LightGlue.git && pip install -e LightGlue
   git clone https://github.com/davnords/LoMa.git && pip install -e LoMa

.. _benchmark_common_options:

Common options
--------------

``--output FILE`` and ``--compare A B``
   ``--output`` saves a run to a JSON file: its command line, its
   environment (the date, the OS, the CPU, the versions of Python, NumPy,
   OpenCV and HomLib, the robust estimators HomLib was built with, and the git
   commit of the scripts), its settings and its results. ``--compare``
   (the ``compare A B`` command for HEB) prints the differences B - A between
   two saved runs: with their standard errors for the shares and the AUCs or
   mAAs (a difference smaller than about twice its standard error could be
   chance), and as ratios for the times and median errors. It prints the commands of the two
   runs and the differences between their environments, and warns if they have
   different settings. The times depend on the host and its load:
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
