ARG UBUNTU_VERSION=22.04
FROM ubuntu:${UBUNTU_VERSION}

ARG COMPILER_VERSION=11

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y -q && \
    apt-get upgrade -y -q && \
    apt-get install -y -q --no-install-recommends \
        ca-certificates \
        cmake \
        gcc-${COMPILER_VERSION} \
        g++-${COMPILER_VERSION} \
        git \
        libeigen3-dev \
        libgtest-dev \
        make \
        software-properties-common && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV CC=/usr/bin/gcc-${COMPILER_VERSION}
ENV CXX=/usr/bin/g++-${COMPILER_VERSION}
