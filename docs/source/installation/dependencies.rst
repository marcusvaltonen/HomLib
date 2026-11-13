.. _dependencies:

Dependencies
============

Python users
------------

To install the Python package you need **Python 3.10 or greater**. Pre-built wheels
are available on PyPI for all major platforms:

.. code-block:: bash

   pip install homlib

No additional dependencies are required — the C++ core is bundled in the wheel.

C++ developers
--------------

To build the core C++ library from source you need:

- **CMake** 3.15 or later
- **Eigen** 3.3 or later — on Ubuntu: ``apt-get install libeigen3-dev``
- **PoseLib** — must be built and installed from `source <https://github.com/PoseLib/PoseLib>`__
- **RansacLib** — included as a git submodule (header-only):

  .. code-block:: bash

     git submodule update --init --recursive

Optional (SupeRANSAC backend)
-----------------------------

To build with the `SupeRANSAC <https://github.com/danini/SuperRANSAC>`__ robust
estimation backend you additionally need:

- **Boost** — on Ubuntu: ``apt-get install libboost-all-dev``
- **OpenCV** 4.x — on Ubuntu: ``apt-get install libopencv-dev``
- **SupeRANSAC** submodule:

  .. code-block:: bash

     git submodule update --init thirdparty/superansac

Then build with ``-DHOMLIB_USE_SUPERANSAC=ON``:

.. code-block:: bash

   CMAKE_ARGS="-DHOMLIB_USE_SUPERANSAC=ON" pip install .

Optional (GC-RANSAC and MAGSAC++ backends)
-------------------------------------------

To build with the GC-RANSAC [Barath2018]_ or MAGSAC++ [Barath2020]_ backends you
need **OpenCV** 4.x. Both use GC-RANSAC's code (MAGSAC++ is built on its
samplers and model classes); MAGSAC++ also needs its own submodule, **gflags**
and **OpenMP**:

.. code-block:: bash

   git submodule update --init thirdparty/graph-cut-ransac thirdparty/magsac

Build with one or both flags:

.. code-block:: bash

   CMAKE_ARGS="-DHOMLIB_USE_GCRANSAC=ON -DHOMLIB_USE_MAGSAC=ON" pip install .

Optional (for examples)
-----------------------

Some gallery examples use additional packages:

- **matplotlib** — plotting
- **opencv-python** — image I/O
- **kornia**, **kornia_moons**, **torch** — feature detection in the real-data example
