FROM fedora:34

# Declaring build variables
ARG CMAKE_VERSION=3.18.6

# Creating in and out directories
RUN mkdir /in
RUN mkdir /out

# System update and packages installation
RUN dnf update -y

# Installing Utilities
RUN dnf install -y wget git make

# CMake installation
RUN wget -O cmake.sh https://github.com/Kitware/CMake/releases/download/v$CMAKE_VERSION/cmake-$CMAKE_VERSION-Linux-x86_64.sh
RUN sh cmake.sh --skip-license --exclude-subdir --prefix=/usr/local/
RUN rm cmake.sh

# Compilers installation
RUN dnf install -y clang

# Setting environment variables
ENV CC=/usr/bin/clang
ENV CXX=/usr/bin/clang++

# No TPLs

# Setting workdir to /in
WORKDIR /in