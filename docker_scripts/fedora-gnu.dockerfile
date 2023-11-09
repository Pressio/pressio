ARG FEDORA_VERSION=latest
FROM fedora:${FEDORA_VERSION}

RUN dnf update -y && \
    dnf install -y \
        cmake \
        g++ \
        gcc \
        git \
        gtest-devel \
        hostname \
        make \
        wget && \
    dnf clean all

ENV CC=/usr/bin/gcc
ENV CXX=/usr/bin/g++
