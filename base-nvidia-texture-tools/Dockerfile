FROM ubuntu:18.04
MAINTAINER Hristina Hristova <hristinaih at gmail.com>

# Install dependencies
RUN apt-get update && \
  apt-get install -y build-essential cmake python python-pip libsm6 libxext6 libxrender-dev imagemagick && \
  pip install opencv-python && \
  rm -rf /var/lib/{apt,dpkg,cache,log}/

RUN apt-get update && \
  apt-get install -y git && \
  rm -rf /var/lib/{apt,dpkg,cache,log}/

RUN cd / && \
  mkdir nvidia-texture-tools && \
  git clone https://github.com/castano/nvidia-texture-tools.git /nvidia-texture-tools/sources && \
  cd /nvidia-texture-tools/sources && \
  git checkout b45560cfc4684fec8a79d812a20780e5d79df9b3 && \
  mkdir build && \
  cd build && \
  cmake ..  && \
  make && \
  make install

ADD ./runner.sh /nvidia-texture-tools/sources
ADD ./3cps_analysis.py /nvidia-texture-tools/sources

WORKDIR /nvidia-texture-tools
VOLUME /nvidia-texture-tools/data
ENTRYPOINT ["/nvidia-texture-tools/sources/runner.sh"]
CMD []
