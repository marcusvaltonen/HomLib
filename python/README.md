# homlib
Python wrappers for C++ core code.

## Install
Use the build script in the project root to build ``libHomLib.a`` which will be linked to the python library.
Then create a virtual environment and install homlib using

```console
pip install .
```

or using ``python -m build``, ``uv build``, etc.

## HPatches tests
Get glue-factory
`uv pip install gluefactory@https://github.com/cvg/glue-factory.git`

`python experiments/hpatches.py`
Having some struggle with `pipeline.run` in that file, not sure why.

## Tests
There are a few examples in the ``example`` folder to demo the wrappers.
