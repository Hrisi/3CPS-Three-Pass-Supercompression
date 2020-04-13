#!/bin/bash

set -e

BASE=/nvidia-texture-tools/data/textures

DATA=${1:-$(cd $BASE && ls *.{jpg,png} || true)}

for f in $DATA; do
  python /nvidia-texture-tools/sources/3cps_analysis.py -i $f -f 0 
done
