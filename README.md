# HomLib

[![Build Status](https://travis-ci.com/marcusvaltonen/HomLib.svg?branch=main)](https://travis-ci.com/marcusvaltonen/HomLib)
![GitHub](https://img.shields.io/github/license/marcusvaltonen/HomLib)
[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/marcusvaltonen/HomLib.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/marcusvaltonen/HomLib/context:cpp)
[![codecov](https://codecov.io/gh/marcusvaltonen/HomLib/branch/main/graph/badge.svg)](https://codecov.io/gh/marcusvaltonen/HomLib)

C++ library for computing homographies with support in MATLAB and Python.

## Solvers
This repository contains four solvers for computing homographies with
simultaneous radial distortion correction and/or incorporating IMU data.

| Authors (year)                | Number of points | Minimal            | Radial distortion coeff. | IMU data           | General homography | Separarate intrinsic/extrinsic |
| ----------------------------- | ---------------- | ------------------ | ------------------------ | ------------------ | ------------------ | ------------------------------ |
| Fitzgibbon (2001)             | 5                |                    | :heavy_check_mark:       |                    | :heavy_check_mark: |                                |
| Kukelova et al. (2015)        | 5                | :heavy_check_mark: | :heavy_check_mark: (2)   |                    | :heavy_check_mark: |                                |
| Valtonen Örnhag et al. (2020) | 4                | :heavy_check_mark: |                          | :heavy_check_mark: |                    | :heavy_check_mark:             |
| Valtonen Örnhag et al. (2021) | 3                | :heavy_check_mark: |                          | :heavy_check_mark: |                    | :heavy_check_mark:             |
|                               | 4                | :heavy_check_mark: | :heavy_check_mark:       | :heavy_check_mark: |                    | :heavy_check_mark:             |

The solvers by Valtonen Örnhag et al. are original implementations, the
other are re-implementations. If you use the code in your work, please cite
the following articles:

```
@InProceedings {valtonen-ornhag-etal-2020-icpr,
    author    = {Valtonen~{\"O}rnhag, Marcus and Persson, Patrik and Wadenb{\"a}ck, M{\aa}rten and {\AA}str{\"o}m, Kalle and Heyden, Anders},
    title     = {Minimal Solvers for Indoor UAV Positioning},
    booktitle = {Proceedings of the 25th International Conference on Pattern Recognition (ICPR)},
    month     = {January},
    year      = {2021},
    pages     = {1136-1143}
}


@InProceedings{valtonen-ornhag-etal-2021-wacv,
    author    = {Valtonen~{\"O}rnhag, Marcus and Persson, Patrik and Wadenb{\"a}ck, M{\aa}rten and {\AA}str{\"o}m, Kalle and Heyden, Anders},
    title     = {Efficient Real-Time Radial Distortion Correction for UAVs},
    booktitle = {Proceedings of the IEEE/CVF Winter Conference on Applications of Computer Vision (WACV)},
    month     = {January},
    year      = {2021},
    pages     = {1751-1760}
}
```

## Dependencies
The source code depends on Eigen 3 (older versions not compatible).
Installation for Ubuntu/Debian:
```bash
    $ apt-get install libeigen3-dev
```
The source code has been compiled and tested on Ubuntu 18.04 (Bionic Beaver) with g++-7 to g++-9 as well
as clang++-7 to clang++-9. Furthermore, it is tested on OSX with Xcode 10-12.

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
