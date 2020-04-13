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
  mkdir nvidia-texture-tool && \
  git clone https://github.com/castano/nvidia-texture-tools.git /nvidia-texture-tool/sources && \
  cd /nvidia-texture-tool/sources && \
  git checkout b45560cfc4684fec8a79d812a20780e5d79df9b3 && \
  mkdir build && \
  cd build && \
  cmake ..  && \
  make && \
  make install

ADD ./runner.sh /nvidia-texture-tool/sources

WORKDIR /nvidia-texture-tool
VOLUME /nvidia-texture-tool/data
ENTRYPOINT ["/nvidia-texture-tool/sources/runner.sh"]
CMD []
