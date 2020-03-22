#!/bin/bash

set -e

BASE=/nvidia-texture-tools/data

DATA=${1:-$(cd $BASE && ls *.{jpg,png} || true)}

for f in $DATA; do
  nvcompress -bc1 $BASE/$f $BASE/${f}.dds 
done
