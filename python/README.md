# homlib
Python wrappers for C++ core code.

## Install
Use the build script in the project root to build ``libHomLib.a`` which will be linked to the python library.
Then create a virtual environment and install homlib using

```console
pip install .
```

or using ``python -m build``, ``uv build``, etc.

## HPatches benchmark
Run the HPatches evaluation with [glue-factory](https://github.com/cvg/glue-factory):

```console
uv pip install gluefactory@https://github.com/cvg/glue-factory.git
python experiments/hpatches.py
```

## Benchmarks
Feature extraction and the benchmarks are separate scripts in ``experiments/``:

```console
uv pip install opencv-python tqdm
python experiments/hpatches_features.py --data_dir data --output data/hpatches_sift.npz
python experiments/hpatches_benchmark.py data/hpatches_sift.npz --output sift.json
python experiments/synthetic_benchmark.py --output synthetic.json
```

- ``hpatches_features.py`` downloads HPatches to ``data`` if needed, and writes the matches of its
  580 pairs: OpenCV SIFT by default, or LoMa with ``--features loma`` (needs ``loma``).
- ``hpatches_benchmark.py`` runs every point solver in LOMSAC on these matches, and reports the share
  of the true inliers found, the AUC of the corner error and the time.
- ``synthetic_benchmark.py`` measures the numerical stability (noise-free), the sensitivity to noise
  and the time of every point solver, on problems from ``generate_problem_instance``.

``--output`` saves the results, and ``--compare before.json after.json`` prints the differences
between two runs (e.g. with two versions of HomLib).

### Test data for the regression test
``tests/test_hpatches.py`` runs the point solvers on the matches of 8 HPatches pairs. The file isn't
in the repository: the test looks for it in ``$HOMLIB_TEST_DATA``, then downloads it from a
release asset, and is skipped without it. To regenerate the file and publish it (maintainers):

```console
python experiments/hpatches_features.py --matches data/hpatches_sift.npz --export_subset hpatches_sift_subset.npz
shasum -a 256 hpatches_sift_subset.npz > hpatches_sift_subset.npz.sha256
gh release create test-data-v1 hpatches_sift_subset.npz hpatches_sift_subset.npz.sha256 \
    --title "Test data" --notes "Data files for the tests (not a HomLib release)"
```

The test checks the downloaded file against the ``.sha256`` file next to it, so a new upload only
needs both files replaced (``gh release upload test-data-v1 ... --clobber``). HPatches is under the
MIT license.

## Tests
There are a few examples in the ``example`` folder to demo the wrappers.
