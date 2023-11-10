ARG FEDORA_VERSION=39
FROM fedora:${FEDORA_VERSION}

RUN dnf update -y && \
    dnf install -y \
        clang \
        cmake \
        git \
        gtest-devel \
        hostname \
        make \
        wget && \
    dnf clean all

ENV CC=/usr/bin/clang
ENV CXX=/usr/bin/clang++
