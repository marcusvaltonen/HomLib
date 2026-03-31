# HomLib
![GitHub release (latest by date)](https://img.shields.io/github/v/release/marcusvaltonen/HomLib)
![GitHub](https://img.shields.io/github/license/marcusvaltonen/HomLib)
![PyPI](https://img.shields.io/pypi/v/homlib)


C++ library for computing homographies with support in MATLAB and Python.

## Solvers
This repository contains the following solvers for computing homographies with
simultaneous radial distortion correction and/or incorporating IMU data.

| Authors (year)                | Number of points | Minimal            | Radial distortion coeff. | IMU data           | General homography | Separarate intrinsic/extrinsic |
| ----------------------------- | ---------------- | ------------------ | ------------------------ | ------------------ | ------------------ | ------------------------------ |
| Fitzgibbon (2001)             | 5                |                    | :heavy_check_mark: (e)   |                    | :heavy_check_mark: |                                |
|                               | 5                |                    | :heavy_check_mark: (1)   |                    | :heavy_check_mark: |                                |
| Kukelova et al. (2015)        | 5                | :heavy_check_mark: | :heavy_check_mark: (2)   |                    | :heavy_check_mark: |                                |
|                               | 6                |                    | :heavy_check_mark: (2)   |                    | :heavy_check_mark: |                                |
| Valtonen Örnhag et al. (2020) | 4                | :heavy_check_mark: |                          | :heavy_check_mark: |                    | :heavy_check_mark:             |
| Valtonen Örnhag et al. (2021) | 3                | :heavy_check_mark: |                          | :heavy_check_mark: |                    | :heavy_check_mark:             |
|                               | 4                | :heavy_check_mark: | :heavy_check_mark: (e)   | :heavy_check_mark: |                    | :heavy_check_mark:             |
| Nakano (2024)                 | 5                | :heavy_check_mark: | :heavy_check_mark: (1)   |                    | :heavy_check_mark: |                                |
| Wadenbäck et al. (2026)       | 5                | :heavy_check_mark: | :heavy_check_mark: (1)   |                    | :heavy_check_mark: |                                |
|                               | 5                | :heavy_check_mark: | :heavy_check_mark: (e)   |                    | :heavy_check_mark: |                                |
|                               | 5                | :heavy_check_mark: | :heavy_check_mark: (2)   |                    | :heavy_check_mark: |                                |

We use the following convention for the different cases: (1) - single-sided, (e) two-sided and equal, and (2) two-sided.

The solvers by Valtonen Örnhag et al. and Wadenbäck et al. are original implementations, the
others are re-implementations. If you use the code in your work, please cite
the respective article:

```
@InProceedings{fitzgibbon-etal-2001-cvpr,
    author    = {Fitzgibbon, Andrew},
    title     = {Simultaneous linear estimation of multiple view geometry and lens distortion},
    booktitle = {Proceedings of the IEEE Computer Society Conference on Computer Vision and Pattern Recognitioni (CVPR)},
    year      = {2001},
}

@InProceedings{kukelova-etal-2015-cvpr,
    author    = {Kukelova, Zuzana and Heller, Jan and Bujnak, Martin and Pajdla, Tomas},
    title     = {Radial Distortion Homography},
    booktitle = {Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition (CVPR)},
    month     = {June},
    year      = {2015}
}

@InProceedings{valtonen-ornhag-etal-2020-icpr,
    author    = {Valtonen~{\"O}rnhag, Marcus and Persson, Patrik and Wadenb{\"a}ck, M{\aa}rten and {\AA}str{\"o}m, Kalle and Heyden, Anders},
    title     = {Minimal Solvers for Indoor {UAV} Positioning},
    booktitle = {Proceedings of the International Conference on Pattern Recognition (ICPR)},
    month     = {January},
    year      = {2021},
    pages     = {1136-1143}
}

@InProceedings{valtonen-ornhag-etal-2021-wacv,
    author    = {Valtonen~{\"O}rnhag, Marcus and Persson, Patrik and Wadenb{\"a}ck, M{\aa}rten and {\AA}str{\"o}m, Kalle and Heyden, Anders},
    title     = {Efficient Real-Time Radial Distortion Correction for {UAV}s},
    booktitle = {Proceedings of the IEEE/CVF Winter Conference on Applications of Computer Vision (WACV)},
    month     = {January},
    year      = {2021},
    pages     = {1751-1760}
}

@InProceedings{nakano-2024-icpr,
    author    = {Nakano, Gaku},
    title     = {Inverse DLT Method for One-Sided Radial Distortion Homography},
    booktitle = {Proceedings of the International Conference on Pattern Recognition (ICPR)},
    year      = {2024},
    pages     = {448-462}
}

@InProceedings{wadenback-etal-2026-3dv,
    author    = {Wadenb{\"a}ck, M{\aa}rten and Valtonen~{\"O}rnhag, Marcus and Edstedt, Johan},
    title     = {Radially Distorted Homographies, Revisited},
    booktitle = {Proceedings of the International Conference on 3D Vision (3DV)},
    year      = {2026},
}
```

## Dependencies
The source code depends on Eigen 3 (older versions not compatible).
Installation for Ubuntu/Debian:
```bash
    $ apt-get install libeigen3-dev
```

Furthermore, [PoseLib](https://github.com/PoseLib/PoseLib/tree/master) is required. Follow the installation
instruction in the repo.

If you want to use the solvers in an LOMSAC framework, we rely
on [RansacLib](https://github.com/tsattler/RansacLib), which is included as a submodule. You can recursively
clone it

```console
git submodule update --init --recursive
```

## Using the solver in MATLAB
(OUTDATED)
It is possible to MEX-compile the solver and use it in MATLAB. Check the
`compile_mex.m` function in the MATLAB directory. You may have to change the path to Eigen,
e.g. `/usr/local/include/eigen3`.

## Using the solver in Python
Note: Solvers with known IMU data is not yet supported in the python package.
```bash
    $ pip install homlib
```

## About the solvers
Many of the solvers were generated using the automatic generator proposed by
Larsson et al. "Efficient Solvers for Minimal Problems by Syzygy-based
Reduction" (CVPR 2017)
