# 3CPS-Three-Pass-Supercompression
This repository contains the code for performing the DXT1 re-compression step of the proposed 3CPS supercompression algorithm.

The re-compression step, based on the Maximum Difference Cut, is implemented in the NVIDIA Texture Tools (repository: https://github.com/castano/nvidia-texture-tools).
A primary feature of NVIDIA Texture Tools is the DXT compression. Our solution contributes to the `nvtt` library, optimizing the `CompressorDXT1` functionaly (including functions for performing DXT1 texture compression).

To install and run our 3CPS solution, you can use the support for building Docker image and running encoding experiments in a Docker container. The folder `textures/` contains uncompressed examplar textures.

## Docker Build
```
$ docker build -t 3cps-nvidia-texture-tools:latest .
```

## Docker Run

To run the software for all examplar textures in `textures/` folder:
```
$ docker run -ti -v $(pwd)/textures:/nvidia-texture-tools/data --rm 3cps-nvidia-texture-tool
```

To run the software for a particular texture `t.jpg` from the `textures/` folder:
```
$ docker run -ti -v $(pwd)/textures:/nvidia-texture-tools/data --rm 3cps-nvidia-texture-tool t.jpg
```

Moreover, to compare the re-compression step from our 3CPS solution to the original DXT1 compression, build the Docker image in the folder `base-nvidia-texture-tools/`. The build will clone the original nvidia-texture-tools repository and install the software. Then, run the Docker container on all examplar textures or on a particular texture.
