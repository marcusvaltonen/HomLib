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

Optional (for examples)
-----------------------

Some gallery examples use additional packages:

- **matplotlib** — plotting
- **opencv-python** — image I/O
- **kornia**, **kornia_moons**, **torch** — feature detection in the real-data example
