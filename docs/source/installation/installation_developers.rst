.. _installation_developers:

Step-by-step installation for developers
========================================

To build the core C++ library you need Eigen 3.3+. On Ubuntu:


.. code-block::

       apt install libeigen3-dev

and also `PoseLib <https://github.com/PoseLib/PoseLib>`__. Follow
the installation instructions in the repo.

Furthermore, you need `RansacLib <https://github.com/tsattler/RansacLib>`__, which
is included as a submodule. You may simply pull it recursively,

.. code-block::

	git submodule update --init --recursive
	
RansacLib is a header-only library and does not need any further installation.

You may now use the **build.sh** script in the root of the repository to build the
C++ library for HomLib.


