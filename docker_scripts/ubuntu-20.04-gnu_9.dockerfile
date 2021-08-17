FROM ubuntu:focal
# Declaring build variables
ARG TZ=Europe/Warsaw
ARG CMAKE_VERSION=3.16.8
ARG GNU_VER=9
ARG GCC=gcc-$GNU_VER
ARG CXX=g++-$GNU_VER

# Setting timezone
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Creating in and out directories
RUN mkdir /in
RUN mkdir /out

# System update and packages installation
RUN apt-get update && apt-get upgrade -y
# Installing Utilities
RUN apt-get install -y wget git make

# CMake installation
RUN wget -O cmake.sh https://github.com/Kitware/CMake/releases/download/v$CMAKE_VERSION/cmake-$CMAKE_VERSION-Linux-x86_64.sh
RUN sh cmake.sh --skip-license --exclude-subdir --prefix=/usr/local/

# Compilers installation
RUN apt-get install -y $GCC $CXX
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/$GCC 10
RUN update-alternatives --install /usr/bin/g++ g++ /usr/bin/$CXX 10
RUN update-alternatives --install /usr/bin/cc cc /usr/bin/gcc 20
RUN update-alternatives --set cc /usr/bin/gcc
RUN update-alternatives --install /usr/bin/c++ c++ /usr/bin/g++ 20
RUN update-alternatives --set c++ /usr/bin/g++

# Setting environment variables
ENV CC=/usr/bin/gcc
ENV CXX=/usr/bin/g++

# No TPLs

# Setting workdir to /in
WORKDIR /in