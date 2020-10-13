# HomLib
C++ library for computing homographies with support in MATLAB and Python.

## Dependencies
The source code depends on Eigen 3 (older versions not compatible).
Installation for Ubuntu/Debian:
```bash
    $ apt-get install libeigen3-dev
```
Tested on Eigen 3.3.4.

## Using the solver in MATLAB
It is possible to MEX-compile the solver and use it in MATLAB. Check the
`compile_mex.m` function in the MATLAB directory. You may have to change the path to Eigen,
e.g. `/usr/local/include/eigen3`.

## Using the solver in Python
This is in the process of being package and will be available on PyPi shortly.

## About the solvers
Many of the solvers were generated using the automatic generator proposed by
Larsson et al. "Efficient Solvers for Minimal Problems by Syzygy-based
Reduction" (CVPR 2017)
