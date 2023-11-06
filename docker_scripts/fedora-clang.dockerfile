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
        make \
        wget && \
    dnf clean all

# Setting environment variables
ENV CC=/usr/bin/clang
ENV CXX=/usr/bin/clang++

# Creating in and out directories
RUN mkdir /in
RUN mkdir /out

# Setting workdir to /in
WORKDIR /in
