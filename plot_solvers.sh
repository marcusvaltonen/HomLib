#!/bin/bash

set -euo pipefail

NOISE=${1:-0}
ITER=${2:-100000}

echo Removing old
rm -f *.csv

./build/example $NOISE $ITER

python tools/plot.py
