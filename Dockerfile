FROM ubuntu:18.04
MAINTAINER Hristina Hristova <hristinaih at gmail.com>

# Install dependencies
RUN apt-get update && \
  apt-get install -y build-essential cmake python python-pip libsm6 libxext6 libxrender-dev imagemagick && \
  pip install opencv-python && \
  rm -rf /var/lib/{apt,dpkg,cache,log}/

ADD ./README.md ./CMakeLists.txt ./LICENSE /nvidia-texture-tools/sources/
ADD ./cmake /nvidia-texture-tools/sources/cmake
ADD ./extern /nvidia-texture-tools/sources/extern
ADD ./src /nvidia-texture-tools/sources/src

RUN cd /nvidia-texture-tools/sources && \
  mkdir build && \
  cd build && \
  cmake ..  && \
  make && \
  make install

ADD ./runner.sh /nvidia-texture-tools/sources
ADD ./base-nvidia-texture/tools/3cps_analysis.py /nvidia-texture-tools/sources

WORKDIR /nvidia-texture-tools
VOLUME /nvidia-texture-tools/data
ENTRYPOINT ["/nvidia-texture-tools/sources/runner.sh"]
CMD []
