FROM fedora:34

# Declaring build variables
ARG CMAKE_VERSION=3.18.6

# Creating in and out directories
RUN mkdir /in
RUN mkdir /out

# System update and packages installation
RUN dnf update -y

# Installing Utilities
RUN dnf install -y wget git make hostname

# CMake installation
RUN wget -O cmake.sh https://github.com/Kitware/CMake/releases/download/v$CMAKE_VERSION/cmake-$CMAKE_VERSION-Linux-x86_64.sh
RUN sh cmake.sh --skip-license --exclude-subdir --prefix=/usr/local/
RUN rm cmake.sh

# Compilers installation
RUN dnf install -y gcc g++

# Setting environment variables
ENV CC=/usr/bin/gcc
ENV CXX=/usr/bin/g++

# Building TPLs
RUN git clone https://github.com/Pressio/pressio-builder.git
WORKDIR /pressio-builder
RUN git checkout main
RUN sed -i 's/source .\/shared\/colors.sh/# colors commnted out/g' main_tpls.sh
RUN chmod +x main_tpls.sh
RUN ./main_tpls.sh -dryrun=no -build-mode=Release -target-dir=/usr/local -tpls=gtest

# Setting workdir to /in
WORKDIR /in