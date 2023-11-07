ARG FEDORA_VERSION=latest
FROM fedora:${FEDORA_VERSION}

RUN dnf update -y && \
    dnf install -y \
        cmake \
        eigen3-devel \
        g++ \
        gcc \
        git \
        gtest-devel \
        hostname \
        make && \
    dnf clean all

ENV CC=/usr/bin/gcc
ENV CXX=/usr/bin/g++
