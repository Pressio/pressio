ARG FEDORA_VERSION=latest
FROM fedora:${FEDORA_VERSION}

RUN dnf update -y && \
    dnf install -y wget git make hostname cmake clang eigen3-devel gtest-devel

# Setting environment variables
ENV CC=/usr/bin/clang
ENV CXX=/usr/bin/clang++

# Creating in and out directories
RUN mkdir /in
RUN mkdir /out

# Setting workdir to /in
WORKDIR /in
