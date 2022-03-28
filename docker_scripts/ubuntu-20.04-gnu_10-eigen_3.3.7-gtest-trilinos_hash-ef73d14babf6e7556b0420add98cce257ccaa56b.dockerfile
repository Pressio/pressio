FROM ubuntu:focal
# Declaring build variables
ARG TZ=Europe/Warsaw
ARG CMAKE_VERSION=3.18.6
ARG GNU_VER=10
ARG CC=gcc-$GNU_VER
ARG CXX=g++-$GNU_VER
ARG GFORTRAN=gfortran-$GNU_VER

# Setting timezone
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# System update and packages installation
RUN apt-get update && apt-get upgrade -y
# Installing Utilities
RUN apt-get install -y wget git make
# Installing OpenMPI
RUN apt-get install -y openmpi-bin openmpi-doc
# Installing BLAS, LAPACK
RUN apt-get install -y libopenblas-dev liblapack-dev

# CMake installation
RUN wget -O cmake.sh https://github.com/Kitware/CMake/releases/download/v$CMAKE_VERSION/cmake-$CMAKE_VERSION-Linux-x86_64.sh
RUN sh cmake.sh --skip-license --exclude-subdir --prefix=/usr/local/
RUN rm cmake.sh

# Compilers installation
RUN apt-get install -y $CC $CXX $GFORTRAN
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/$CC 10
RUN update-alternatives --install /usr/bin/g++ g++ /usr/bin/$CXX 10
RUN update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/$GFORTRAN 10
RUN update-alternatives --install /usr/bin/cc cc /usr/bin/gcc 20
RUN update-alternatives --set cc /usr/bin/gcc
RUN update-alternatives --install /usr/bin/c++ c++ /usr/bin/g++ 20
RUN update-alternatives --set c++ /usr/bin/g++
RUN update-alternatives --install /usr/bin/fortrann fortrann /usr/bin/gfortran 20
RUN update-alternatives --set fortrann /usr/bin/gfortran

# Setting environment variables
ENV CC=/usr/bin/mpicc
ENV CXX=/usr/bin/mpic++
ENV FC=/usr/bin/mpifort
ENV F77=/usr/bin/mpifort
ENV F90=/usr/bin/mpifort
ENV MPIRUNe=/usr/bin/mpirun

# Building TPLs
WORKDIR /home
RUN mkdir pressio_builds
RUN mkdir pressio_repos
WORKDIR /home/pressio_repos
RUN git clone https://github.com/Pressio/pressio-builder.git
WORKDIR /home/pressio_repos/pressio-builder
RUN git checkout main
RUN sed -i 's/source .\/shared\/colors.sh/# colors commnted out/g' main_tpls.sh
RUN sed -i 's/9fec35276d846a667bc668ff4cbdfd8be0dfea08/ef73d14babf6e7556b0420add98cce257ccaa56b/g' tpls/tpls_versions_details.sh
RUN sed -i 's/GTESTBRANCH=master/GTESTBRANCH=main/g' tpls/tpls_versions_details.sh
RUN chmod +x main_tpls.sh
RUN ./main_tpls.sh -dryrun=no -build-mode=Release -target-dir=../../pressio_builds -tpls=gtest,eigen,trilinos -cmake-generator-names=default,default,default
RUN git reset --hard origin/main

# Cleaning after builds
WORKDIR /home
RUN rm -rf pressio_builds/gtest/build && rm -rf pressio_builds/gtest/googletest
RUN rm -rf pressio_builds/trilinos/build && rm -rf pressio_builds/trilinos/Trilinos
RUN rm -rf pressio_builds/eigen/eigen-3.3.7 && rm pressio_builds/eigen/eigen-3.3.7.tar.gz

# Setting workdir to /
WORKDIR /