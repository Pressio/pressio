FROM ubuntu:focal
# Declaring build variables
ARG TZ=Europe/Warsaw
ARG CMAKE_VERSION=3.18.6

# Setting timezone
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Creating in and out directories
RUN mkdir /in
RUN mkdir /out

# System update and packages installation
RUN apt-get update && apt-get upgrade -y
# Installing Utilities
RUN apt-get install -y wget git make gnupg2 software-properties-common

# CMake installation
RUN wget -O cmake.sh https://github.com/Kitware/CMake/releases/download/v$CMAKE_VERSION/cmake-$CMAKE_VERSION-Linux-x86_64.sh
RUN sh cmake.sh --skip-license --exclude-subdir --prefix=/usr/local/
RUN rm cmake.sh

# Compilers installation
RUN wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN add-apt-repository "deb https://apt.repos.intel.com/oneapi all main"
RUN apt-get install -y intel-hpckit

# Setting environment variables
ENV CC=/opt/intel/oneapi/compiler/latest/linux/bin/intel64/icc
ENV CXX=/opt/intel/oneapi/compiler/latest/linux/bin/intel64/icpc

# No TPLs

# Setting workdir to /in
WORKDIR /in