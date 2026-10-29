# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

HomLib is a C++ library for computing homographies with simultaneous radial distortion correction, optionally incorporating IMU data. It supports point-based, affine-covariant, and orientation-covariant features. Python bindings are provided via pybind11 and published to PyPI as `homlib`.

## Build Commands

### C++ (full build + tests)
```bash
./build.sh
```
Or manually:
```bash
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --config Release -- -j $(nproc)
cd tests && ctest -j $(nproc) --output-on-failure
```

### Run a single C++ test
Tests use Catch2 v3 and are registered via `catch_discover_tests`. After building:
```bash
cd build/tests
ctest -R "test_name_pattern" --output-on-failure
```

### Python virtual environment
Always use the project venv at `venv/` in the project root:
```bash
source venv/bin/activate
```

### Python package (local dev build)
Requires PoseLib and Eigen3 installed. Uses scikit-build-core. Prefer `uv` for faster rebuilds:
```bash
source venv/bin/activate
CMAKE_PREFIX_PATH=/home/elramav/Scripts/PoseLib/install uv pip install .
```

### Python tests
```bash
source venv/bin/activate
pytest -s python/tests
```

## Dependencies

- **Eigen3** (required): `apt-get install libeigen3-dev`
- **PoseLib** (required): must be built and installed from source (see PoseLib GitHub repo). CI pins commit `a69263d5d8`.
- **RansacLib** (submodule, for LOMSAC/RANSAC): `git submodule update --init --recursive`
- **Catch2 v3** (fetched automatically by CMake for tests)
- **pybind11** (fetched automatically by CMake for Python bindings)

## Architecture

### Core library (`HomLib/`)
- **`HomLib/solvers/`** - Each solver lives in a subdirectory named `{author}_{venue}_{year}/`. Entry points are `get_*()` functions returning `std::vector<HomLib::PoseData>`. Some solvers have internal `solver_*` files with elimination template machinery.
- **`HomLib/includes/HomLib/`** - Public headers. Each solver header (e.g. `get_fitzgibbon_cvpr_2001.hpp`) declares the free functions and also CRTP-based `PoseEstimator` subclasses used for RANSAC integration.
- **`HomLib/helpers/`** - Shared utilities: 2D point normalization, Gauss-Jordan elimination, polynomial root finding, radial distortion/undistortion, problem instance generation.
- **`HomLib/robust/`** - RANSAC estimator adapters (`RansacEstimator`, `AffineRansacEstimator`, `OrientationRansacEstimator`) and nonlinear refinement functions.
- **`posedata.hpp`** - Central data struct: `{homography, focal_length, distortion_parameter, distortion_parameter2}`.

### Python bindings (`python/`)
- `python/src/main.cpp` - Single pybind11 module `_core` wrapping all solvers (both direct `estimate_*` and LOMSAC `lomsac_*` variants). Also exposes `PoseData`, `LORansacOptions`, `RansacStatistics`, `ProblemConfig`, `ProblemInstance`, and `generate_problem_instance`.
- `python/src/homlib/__init__.py` - Re-exports everything from `_core`.
- `python/tests/` - Uses pytest with approval tests (expected outputs in `approved_files/`).

### Naming conventions
Solver functions and files follow the pattern `{author}_{venue}_{year}`. Distortion variants use suffixes: `_one_sided`, `_two_sided_equal`, `_two_sided`, `_6pt`. The Python API mirrors this exactly.

### CI
GitHub Actions (`.github/workflows/main.yml`) builds and tests on Linux, macOS, and Windows. On tag push (`v*`), it builds wheels via cibuildwheel and publishes to PyPI.

## Adding a New Solver

1. Create `HomLib/solvers/{author}_{venue}_{year}/` with implementation `.cpp` files
2. Add a public header in `HomLib/includes/HomLib/` declaring `get_*()` functions and a `PoseEstimator` subclass if RANSAC support is needed
3. Register source files in `HomLib/CMakeLists.txt` and add include dirs in the root `CMakeLists.txt`
4. Add C++ tests in `tests/` (Catch2 TEST_CASE with known input/output pairs)
5. Add pybind11 wrappers in `python/src/main.cpp` and export in `python/src/homlib/__init__.py`
