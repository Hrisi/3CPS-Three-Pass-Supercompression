# 3CPS-Three-Pass-Supercompression
This repository contains the code for performing the DXT1 re-compression of the proposed 3CPS supercompression

Work in progress

The re-compression step, based on the Maximum Difference Cut, is implemented in the NVIDIA Texture Tools (repository: https://github.com/castano/nvidia-texture-tools).
A primary feature of NVIDIA Texture Tools is the DXT compression. Our contribution refers to nvtt library, and more specifically, CompressorDXT1 file, containing the functions for performing DXT1 texture compression.

How to build (Linux/OSX):
Use cmake and the provided configure script:

$./configure
$make
$sudo make install

How to run it:
nvcompress -bc1 'out'.png 'in'.dss
