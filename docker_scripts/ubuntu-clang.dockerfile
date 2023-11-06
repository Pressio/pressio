ARG UBUNTU_VERSION=latest
FROM ubuntu:${UBUNTU_VERSION}

ENV CC=/usr/bin/clang
ENV CXX=/usr/bin/clang++
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y -q && \
    apt-get upgrade -y -q && \
    apt-get install -y -q --no-install-recommends \
        ca-certificates \
        clang \
        cmake \
        git \
        libeigen3-dev \
        libgtest-dev \
        make \
        wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Creating in and out directories
RUN mkdir /in
RUN mkdir /out

# Setting workdir to /in
WORKDIR /in
