# HomLib

[![Build Status](https://travis-ci.com/marcusvaltonen/HomLib.svg?branch=main)](https://travis-ci.com/marcusvaltonen/HomLib)
![GitHub](https://img.shields.io/github/license/marcusvaltonen/HomLib)

C++ library for computing homographies with support in MATLAB and Python.
More solvers and documentation will be added soon.

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
The official python repository is [python-homlib](https://github.com/marcusvaltonen/python-homlib).
A pre-alpha release is available at PyPi, and can be installed using
```bash
    $ pip install homlib
```

## About the solvers
Many of the solvers were generated using the automatic generator proposed by
Larsson et al. "Efficient Solvers for Minimal Problems by Syzygy-based
Reduction" (CVPR 2017)
