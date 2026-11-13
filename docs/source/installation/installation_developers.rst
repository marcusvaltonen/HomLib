.. _installation_developers:

Step-by-step installation for developers
========================================

C++ library
-----------

To build the core C++ library you need Eigen 3.3+. On Ubuntu:

.. code-block:: bash

   apt install libeigen3-dev

and also `PoseLib <https://github.com/PoseLib/PoseLib>`__. Follow
the installation instructions in the repo.

Furthermore, you need `RansacLib <https://github.com/tsattler/RansacLib>`__, which
is included as a submodule. You may simply pull it recursively:

.. code-block:: bash

   git submodule update --init --recursive

RansacLib is a header-only library and does not need any further installation.

You may now use the **build.sh** script in the root of the repository to build the
C++ library for HomLib:

.. code-block:: bash

   ./build.sh

Or manually:

.. code-block:: bash

   mkdir -p build && cd build
   cmake -DCMAKE_BUILD_TYPE=Release ..
   cmake --build . --config Release -- -j $(nproc)

Building with SupeRANSAC (optional)
------------------------------------

To enable the `SupeRANSAC <https://github.com/danini/SuperRANSAC>`__ robust
estimation backend, you additionally need Boost and OpenCV 4.x:

.. code-block:: bash

   apt install libboost-all-dev libopencv-dev

Pull the SupeRANSAC submodule:

.. code-block:: bash

   git submodule update --init thirdparty/superansac

Then build with ``-DHOMLIB_USE_SUPERANSAC=ON``:

.. code-block:: bash

   mkdir -p build && cd build
   cmake -DCMAKE_BUILD_TYPE=Release -DHOMLIB_USE_SUPERANSAC=ON ..
   cmake --build . --config Release -- -j $(nproc)

Building with GC-RANSAC and MAGSAC++ (optional)
-------------------------------------------------

To enable the GC-RANSAC [Barath2018]_ or MAGSAC++ [Barath2020]_ backends you
need OpenCV 4.x:

.. code-block:: bash

   apt install libopencv-dev

Pull the graph-cut-ransac submodule, used by both backends (MAGSAC++ is built
on GC-RANSAC's samplers and model classes), and the magsac submodule for
MAGSAC++, which also needs gflags and OpenMP:

.. code-block:: bash

   git submodule update --init thirdparty/graph-cut-ransac thirdparty/magsac

Then build with ``-DHOMLIB_USE_GCRANSAC=ON`` and/or ``-DHOMLIB_USE_MAGSAC=ON``
(each one works without the other):

.. code-block:: bash

   mkdir -p build && cd build
   cmake -DCMAKE_BUILD_TYPE=Release \
         -DHOMLIB_USE_GCRANSAC=ON \
         -DHOMLIB_USE_MAGSAC=ON ..
   cmake --build . --config Release -- -j $(nproc)

All optional backends can be combined in a single build:

.. code-block:: bash

   cmake -DCMAKE_BUILD_TYPE=Release \
         -DHOMLIB_USE_SUPERANSAC=ON \
         -DHOMLIB_USE_GCRANSAC=ON \
         -DHOMLIB_USE_MAGSAC=ON ..

Python package from source
--------------------------

Install the Python package in your virtual environment:

.. code-block:: bash

   pip install .

With optional backends:

.. code-block:: bash

   # SupeRANSAC
   CMAKE_ARGS="-DHOMLIB_USE_SUPERANSAC=ON" pip install .

   # GC-RANSAC and MAGSAC++
   CMAKE_ARGS="-DHOMLIB_USE_GCRANSAC=ON -DHOMLIB_USE_MAGSAC=ON" pip install .

   # All backends
   CMAKE_ARGS="-DHOMLIB_USE_SUPERANSAC=ON -DHOMLIB_USE_GCRANSAC=ON -DHOMLIB_USE_MAGSAC=ON" pip install .

Running tests
-------------

**C++ tests** (Catch2):

.. code-block:: bash

   cd build/tests
   ctest --output-on-failure

Run a single test by name:

.. code-block:: bash

   ctest -R "Wadenback" --output-on-failure

**Python tests** (pytest):

.. code-block:: bash

   pip install pytest pytest-approvaltests numpy
   pytest -s python/tests

Building documentation locally
------------------------------

.. code-block:: bash

   cd docs
   pip install -r requirements.txt
   make html

The built documentation will be in ``docs/_build/html/``.
