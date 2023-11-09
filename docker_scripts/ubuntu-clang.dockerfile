ARG UBUNTU_VERSION=22.04
FROM ubuntu:${UBUNTU_VERSION}

ARG COMPILER_VERSION=14

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y -q && \
    apt-get upgrade -y -q && \
    apt-get install -y -q --no-install-recommends \
        ca-certificates \
        clang-${COMPILER_VERSION} \
        cmake \
        git \
        libeigen3-dev \
        libgtest-dev \
        make && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV CC=/usr/bin/clang-${COMPILER_VERSION}
ENV CXX=/usr/bin/clang++-${COMPILER_VERSION}
