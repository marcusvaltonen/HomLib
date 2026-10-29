# HomLib
![GitHub release (latest by date)](https://img.shields.io/github/v/release/marcusvaltonen/HomLib)
![GitHub](https://img.shields.io/github/license/marcusvaltonen/HomLib)
![PyPI](https://img.shields.io/pypi/v/homlib)
[![Documentation Status](https://readthedocs.org/projects/homlib/badge/?version=latest)](https://homlib.readthedocs.io/en/latest/?badge=latest)


C++ library for computing homographies with support in Python.

## Quick start

```python
import numpy as np
import homlib

# Estimate homography with one-sided radial distortion
solver = homlib.Wadenback3DV2026OneSided()
poses = solver.estimate(x, y)  # returns list[PoseData]

# Robust estimation with LOMSAC
lomsac = homlib.LOMSAC(solver, squared_inlier_threshold=0.005)
estimate, stats = lomsac.estimate(x_noisy, y_noisy)
```

If built with SupeRANSAC support:

```python
sa = homlib.SupeRANSAC(solver, inlier_threshold=0.07)
estimate, stats = sa.estimate(x_noisy, y_noisy)
```

## Solvers
This repository contains the following solvers for computing homographies with
simultaneous radial distortion correction and/or incorporating IMU data.

We further support affine-covariant and orientation-covariant solvers, requiring
fewer correspondences than their point-based counterparts by leveraging the geometric
information of the corresponding descriptors.

| Authors (year)                           | #Corrs | Minimal            | Radial distortion coeff. | IMU data           | General homography | Descriptor |
| ---------------------------------------- | ------ | ------------------ | ------------------------ | ------------------ | ------------------ | ---------- |
| Fitzgibbon (2001)                        | 5      |                    | :heavy_check_mark: (e)   |                    | :heavy_check_mark: | P          |
|                                          | 5      |                    | :heavy_check_mark: (1)   |                    | :heavy_check_mark: | P          |
| Kukelova et al. (2015)                   | 5      | :heavy_check_mark: | :heavy_check_mark: (2)   |                    | :heavy_check_mark: | P          |
|                                          | 6      |                    | :heavy_check_mark: (2)   |                    | :heavy_check_mark: | P          |
| Valtonen Örnhag et al. (2020)            | 4      | :heavy_check_mark: |                          | :heavy_check_mark: |                    | P          |
| Valtonen Örnhag et al. (2021)            | 3      | :heavy_check_mark: |                          | :heavy_check_mark: |                    | P          |
|                                          | 4      | :heavy_check_mark: | :heavy_check_mark: (e)   | :heavy_check_mark: |                    | P          |
| Nakano (2024)                            | 5      | :heavy_check_mark: | :heavy_check_mark: (1)   |                    | :heavy_check_mark: | P          |
| Wadenbäck et al. (2026)                  | 5      | :heavy_check_mark: | :heavy_check_mark: (1)   |                    | :heavy_check_mark: | P          |
|                                          | 5      | :heavy_check_mark: | :heavy_check_mark: (e)   |                    | :heavy_check_mark: | P          |
|                                          | 5      | :heavy_check_mark: | :heavy_check_mark: (2)   |                    | :heavy_check_mark: | P          |
| Barath and Hajder (2016)                 | 2      |                    |                          |                    | :heavy_check_mark: | A          |
| Valtonen Örnhag and Adalbjörnsson (2026) | 2      |                    | :heavy_check_mark: (1)   |                    | :heavy_check_mark: | A          |
|                                          | 3      | :heavy_check_mark: | :heavy_check_mark: (1)   |                    | :heavy_check_mark: | O          |

We use the following convention for the different distortion cases: (1) - single-sided, (e) two-sided and equal, and (2) two-sided.
The following descriptor types: (P) - point-based, (A) affine-covariant, and (O) orientation-covariant.

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

@InProceedings{valtonen-ornhag-adalbjornsson-2026-icpr,
    author    = {Valtonen~{\"O}rnhag, Marcus and Adalbj{\"o}rnsson, Stefan},
    title     = {Radial Distortion Homography Estimation from Affine-Covariant or Orientation-Covariant Features},
    booktitle = {Proceedings of the International Conference on Pattern Recognition (ICPR)},
    year      = {2026},
}

@Article{barath-2025-superansac,
    author    = {Barath, Daniel},
    title     = {{SupeRANSAC}: One {RANSAC} to Rule Them All},
    journal   = {arXiv preprint arXiv:2506.04803},
    year      = {2025},
}
```

## Installation

### Python (from PyPI)
Pre-built wheels for all standard platforms and supported Python versions:
```bash
pip install homlib
```

### Python (from source)
Requires [PoseLib](https://github.com/PoseLib/PoseLib/tree/master) and Eigen 3.

```bash
pip install .
```

#### With SupeRANSAC support
Requires Boost and OpenCV in addition to the above.

```bash
git submodule update --init thirdparty/superansac
CMAKE_ARGS="-DHOMLIB_USE_SUPERANSAC=ON" pip install .
```

### C++ (from source)
Dependencies:
- **Eigen 3**: `apt-get install libeigen3-dev`
- **[PoseLib](https://github.com/PoseLib/PoseLib)**: build and install from source
- **RansacLib** (submodule): `git submodule update --init --recursive`

```bash
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --config Release
```

## Python usage

### Direct estimation
```python
import homlib

# Point-based solver
solver = homlib.Wadenback3DV2026OneSided()
poses = solver.estimate(x, y)  # x, y are 2xN numpy arrays

# Affine-covariant solver (fewer correspondences needed)
aff_solver = homlib.ValtonenOrnhagICPR2026AffineSingleSided()
poses = aff_solver.estimate(x, y, A=affine_matrices)  # A is 4xN

# Orientation-covariant solver
ori_solver = homlib.ValtonenOrnhagICPR2026OrientationSingleSided()
poses = ori_solver.estimate(x, y, ori=orientations)  # ori is 2xN
```

### Robust estimation with LOMSAC
```python
solver = homlib.Wadenback3DV2026OneSided()
lomsac = homlib.LOMSAC(solver, squared_inlier_threshold=0.005)
estimate, stats = lomsac.estimate(x, y)

print(estimate.homography)
print(estimate.distortion_parameter)
print(f"{len(stats.inlier_indices)} inliers")
```

### Robust estimation with SupeRANSAC
Requires building with `-DHOMLIB_USE_SUPERANSAC=ON`.

```python
solver = homlib.Wadenback3DV2026OneSided()
sa = homlib.SupeRANSAC(solver, inlier_threshold=0.07)
estimate, stats = sa.estimate(x, y)

print(f"{len(stats.inliers)} inliers, {stats.num_iterations} iterations")
```

Any solver can be used with any RANSAC backend:
```python
for solver in [homlib.FitzgibbonCVPR2001OneSided(),
               homlib.NakanoICPR2025OneSided(),
               homlib.Wadenback3DV2026OneSided()]:
    lomsac = homlib.LOMSAC(solver, squared_inlier_threshold=0.005)
    estimate, stats = lomsac.estimate(x, y)
```

## About the solvers
Many of the solvers were generated using the automatic generator proposed by
Larsson et al. "Efficient Solvers for Minimal Problems by Syzygy-based
Reduction" (CVPR 2017)
