FROM ubuntu:focal
# Declaring build variables
ARG TZ=Europe/Warsaw
ARG CMAKE_VERSION=3.18.6
ARG CLANG_VER=9
ARG CC=clang-$CLANG_VER
ARG CXX=clang++-$CLANG_VER

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
RUN rm cmake.sh

# Compilers installation
# NOTE: clang-12 installs libstdc++-9-dev on focal
#       ending up without c++20 headers like <concepts>
#       so request libstdc++-10-dev explicitly
RUN apt-get install -y $CC $CXX libstdc++-10-dev

# Setting environment variables
ENV CC=/usr/bin/clang-$CLANG_VER
ENV CXX=/usr/bin/clang++-$CLANG_VER

# Building TPLs
RUN git clone https://github.com/Pressio/pressio-builder.git
WORKDIR /pressio-builder
RUN git checkout main
RUN sed -i 's/source .\/shared\/colors.sh/# colors commnted out/g' main_tpls.sh
RUN chmod +x main_tpls.sh
RUN ./main_tpls.sh -dryrun=no -build-mode=Release -target-dir=/usr/local -tpls=eigen,gtest

# Setting workdir to /in
WORKDIR /in