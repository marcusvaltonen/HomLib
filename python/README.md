# homlib
Python wrappers for C++ core code.

## Install
Use the build script in the project root to build ``libHomLib.a`` which will be linked to the python library.
Then create a virtual environment and install homlib using

```console
pip install .
```

or using ``python -m build``, ``uv build``, etc.

## Benchmarks
The scripts in ``experiments/`` compare the solvers on real image pairs (HPatches) and on synthetic
problems. Computing the matches of the image pairs is a separate step, so that every later run
compares the solvers on the same matches:

```console
uv pip install opencv-python tqdm
python experiments/hpatches_features.py --data_dir data --output data/hpatches_sift.npz
python experiments/hpatches_benchmark.py data/hpatches_sift.npz --output sift.json
python experiments/synthetic_benchmark.py --output synthetic.json
```

- ``hpatches_features.py`` downloads HPatches to ``data`` if needed (a large archive), and writes the
  matches of its 580 pairs. It uses OpenCV SIFT by default, with precise upscaling (OpenCV 4.8 or
  later) and Lowe's ratio test. It can also use LightGlue with SuperPoint, DISK, ALIKED or SIFT
  (``--matcher lightglue``), or LoMa (``--features loma``). Wadenbäck et al. (3DV 2026) use SuperPoint with LightGlue on HPatches,
  and the docstring of ``hpatches_features.py`` gives their settings. LightGlue and LoMa install
  from https://github.com/cvg/LightGlue and https://github.com/davnords/LoMa; the PyPI package
  named ``loma`` is unrelated.
- ``hpatches_benchmark.py`` runs every point solver in LOMSAC on these matches. It reports the share
  of the true inliers found, the AUC of the corner error, and the time.
- ``synthetic_benchmark.py`` runs every point solver on problems from ``generate_problem_instance``.
  It reports how often each solver misses the true solution, how its error grows with noise, and
  its time. With ``--outliers``, it runs LOMSAC on problems with outliers instead, and reports the share of
  the true inliers that LOMSAC finds.

``--output`` saves the results, and ``--compare before.json after.json`` prints the differences
between two runs, for example with two versions of HomLib.

``experiments/hpatches.py`` runs another HPatches evaluation, with
[glue-factory](https://github.com/cvg/glue-factory):

```console
uv pip install gluefactory@https://github.com/cvg/glue-factory.git
python experiments/hpatches.py
```

### Test data for the regression test
``tests/test_hpatches.py`` runs the point solvers on the SIFT matches of 8 HPatches pairs. These
matches aren't in the repository: the test reads them from ``$HOMLIB_TEST_DATA``, or downloads them
from a release asset (``$HOMLIB_TEST_DATA_URL``, for a fork). Without them, the test is skipped
locally, and fails in CI.

So the release must exist before CI runs this test. To publish it (maintainers):

```console
python experiments/hpatches_features.py --matches data/hpatches_sift.npz --export_subset hpatches_sift_subset.npz
shasum -a 256 hpatches_sift_subset.npz > hpatches_sift_subset.npz.sha256
gh release create test-data-v1 hpatches_sift_subset.npz hpatches_sift_subset.npz.sha256 \
    -R marcusvaltonen/HomLib --latest=false \
    --title "Test data" --notes "Data files for the tests (not a HomLib release)"
```

To change the data, publish it under a new tag (``test-data-v2``) and update ``SUBSET_URL`` in the
test: existing checkouts and caches keep the old file. HPatches is under the MIT license, which
allows publishing these matches.

## Tests
There are a few examples in the ``example`` folder to demo the wrappers.
