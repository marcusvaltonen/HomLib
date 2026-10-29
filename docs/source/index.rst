Welcome to homlib's documentation!
=========================================================

**HomLib** is a C++ library (with Python bindings) for computing homographies with
simultaneous radial distortion correction, optionally incorporating IMU data.
It supports point-based, affine-covariant, and orientation-covariant features.

.. code-block:: python

   import numpy as np
   import homlib

   # Generate a minimal 5-point sample
   H_gt = np.random.randn(3, 3)
   x = np.random.randn(2, 5)
   y = H_gt @ np.vstack((x, np.ones((1, 5))))
   y = y[:2] / y[2]

   # Estimate homography with one-sided radial distortion
   poses = homlib.estimate_wadenback_3dv_2026_one_sided(x, y, True)
   print(poses[0].homography)

Install from PyPI:

.. code-block:: bash

   pip install homlib

| `GitHub <https://github.com/marcusvaltonen/HomLib>`__ · `PyPI <https://pypi.org/project/homlib/>`__ · :ref:`api`

**Maintainer**: `Marcus Valtonen Örnhag <marcusvaltonen@gmail.com>`_

Citing HomLib
-------------

If you use HomLib in your research, please cite the relevant paper(s) from
the :ref:`references` page. For the library itself:

.. code-block:: bibtex

   @misc{HomLib,
       title  = {HomLib -- Homography Library of state-of-the-art methods},
       author = {Marcus Valtonen {\"O}rnhag},
       URL    = {https://github.com/marcusvaltonen/HomLib},
       year   = {2020},
   }


.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Getting started:

   installation/dependencies.rst
   installation/installation_users.rst
   installation/installation_developers.rst

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Examples and usage:

   gallery/index.rst

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Reference documentation:

   api/index.rst
   references.rst
