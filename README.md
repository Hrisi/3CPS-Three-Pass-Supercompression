# 3CPS-Three-Pass-Supercompression
This repository contains the code related to the 3CPS supercompression proposal, which is introduced in
Hristina Hristova, Gwendal Simon, Stefano Petrangeli, and Vishy Swaminathan. ''**3CPS: A Novel Supercompression for the Delivery of 3D Object Textures**, in ACM Multimedia Systems (MMSys) Conference, 2020.

It performs a three-pass compression based on DXT1. The re-compression step, based on the Maximum Difference Cut, is implemented in the NVIDIA Texture Tools (repository: https://github.com/castano/nvidia-texture-tools).
A primary feature of NVIDIA Texture Tools is the DXT compression. Our solution contributes to the `nvtt` library, optimizing the `CompressorDXT1` functionality (including functions for performing DXT1 texture compression).

To install and run our 3CPS solution, you can use the support for building Docker image and running encoding experiments in a Docker container. The folder `experimental_data/textures/` contains uncompressed examplar textures.

## Docker Build
```
$ docker build -t 3cps-nvidia-texture-tools:latest .
```

## Docker Run
Always run the software from the root folder of the `git` repository.

To run the software for all examplar textures in `experimental_data/textures/` folder:
```
$ docker run -ti -v $(pwd)/experimental_data:/nvidia-texture-tools/data --rm 3cps-nvidia-texture-tools
```

To run the software for a particular texture `t.jpg` from the `experimental_data/textures/` folder:
```
$ docker run -ti -v $(pwd)/experimental_data:/nvidia-texture-tools/data --rm 3cps-nvidia-texture-tools t.jpg
```

Moreover, to compare the re-compression step from our 3CPS solution to the original DXT1 compression, build the Docker image in the folder `base-nvidia-texture-tools/`. The build will clone the original nvidia-texture-tools repository and install the software. Then, run the Docker container for all examplar textures or a particular texture.

The compressed images are stored in folders `experimental_data/3cps/` (for the 3CPS supercompression method) and `experimental_data/dxt1/` (for the DXT1 method from the NVIDIA Texture Tools).
