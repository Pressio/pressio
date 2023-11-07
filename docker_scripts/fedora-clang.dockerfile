ARG FEDORA_VERSION=latest
FROM fedora:${FEDORA_VERSION}

RUN dnf update -y && \
    dnf install -y \
        clang \
        cmake \
        eigen3-devel \
        git \
        gtest-devel \
        hostname \
        make && \
    dnf clean all

ENV CC=/usr/bin/clang
ENV CXX=/usr/bin/clang++
