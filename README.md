# 3CPS-Three-Pass-Supercompression
This repository contains the code for performing the DXT1 re-compression of the proposed 3CPS supercompression

The re-compression step, based on the Maximum Difference Cut, is implemented in the NVIDIA Texture Tools (repository: https://github.com/castano/nvidia-texture-tools).
A primary feature of NVIDIA Texture Tools is the DXT compression. Our contribution refers to nvtt library, and more specifically, CompressorDXT1 file, containing the functions for performing DXT1 texture compression.

There is support for building Docker image and running the encoding experiments in a Docker container.

## Docker Build
```
$ docker build -t 3cps-nvidia-texture-tools:latest .
```

## Docker Run

To run the software for all examplar textures in `textures/` folder:
```
$ docker run -ti -v $(pwd)/textures:/nvidia-texture-tools/data --rm 3cps-nvidia-texture-tool
```

To run the software for all particular texture `t.jpg` from the `textures/` folder:
```
$ docker run -ti -v $(pwd)/textures:/nvidia-texture-tools/data --rm 3cps-nvidia-texture-tool t.jpg
```
