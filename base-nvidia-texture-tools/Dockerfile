FROM ubuntu:18.04
MAINTAINER Hristina Hristova <hristinaih at gmail.com>

# Install dependencies
RUN apt-get update && \
  apt-get install -y build-essential cmake && \
  rm -rf /var/lib/{apt,dpkg,cache,log}/

RUN apt-get update && \
  apt-get install -y git && \
  rm -rf /var/lib/{apt,dpkg,cache,log}/

RUN cd / && \
  mkdir nvidia-texture-tools && \
  git clone https://github.com/castano/nvidia-texture-tools.git /nvidia-texture-tools/sources && \
  cd /nvidia-texture-tools/sources && \
  mkdir build && \
  cd build && \
  cmake ..  && \
  make && \
  make install

ADD ./runner.sh /nvidia-texture-tools/sources

WORKDIR /nvidia-texture-tools
VOLUME /nvidia-texture-tools/data
ENTRYPOINT ["/nvidia-texture-tools/sources/runner.sh"]
CMD []
