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
The scripts in ``experiments/`` compare the solvers and robust estimators on real image pairs
(HPatches) and on synthetic problems, and ``tests/test_hpatches.py`` is a regression test on 8
HPatches pairs:

```console
uv pip install opencv-python tqdm
python experiments/hpatches_features.py --data_dir data --output data/hpatches_sift.npz
python experiments/hpatches_benchmark.py data/hpatches_sift.npz --output sift.json
python experiments/synthetic_benchmark.py --output synthetic.json
```

The [Benchmarking guide](https://homlib.readthedocs.io/en/latest/guides/benchmarking.html)
([source](../docs/source/guides/benchmarking.rst)) explains what they measure, their options, how
to compare two runs, and how to publish the test data of the regression test.

## Tests
There are a few examples in the ``example`` folder to demo the wrappers.
