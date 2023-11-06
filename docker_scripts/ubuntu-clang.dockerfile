ARG UBUNTU_VERSION=latest
FROM ubuntu:${UBUNTU_VERSION}
# Declaring build variables
ARG CMAKE_VERSION=3.18.6

# Creating in and out directories
RUN mkdir /in
RUN mkdir /out

ENV CC=/usr/bin/clang
ENV CXX=/usr/bin/clang++
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y -q && \
    apt-get upgrade -y -q && \
    apt-get install -y -q --no-install-recommends \
        ca-certificates \
        clang \
        git \
        libeigen3-dev \
        libgtest-dev \
        make \
        wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# CMake installation
RUN wget -O cmake.sh https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-Linux-x86_64.sh
RUN sh cmake.sh --skip-license --exclude-subdir --prefix=/usr/local/
RUN rm cmake.sh

# Setting workdir to /in
WORKDIR /in
