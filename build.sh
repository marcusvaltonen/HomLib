#!/bin/bash

set -euo pipefail

mkdir -p build && cd build

# Configure
cmake -DCMAKE_BUILD_TYPE=Release ..

# Build (for Make on Unix equivalent to `make -j $(nproc)`)
cmake --build . --config Release -- -j $(nproc)

# Test
cd tests
ctest -j $(nproc) --output-on-failure
cd ..
