FROM fedora:34

# Declaring build variables
ARG CMAKE_VERSION=3.16.8

# Creating in and out directories
RUN mkdir /in
RUN mkdir /out

# System update and packages installation
RUN dnf update -y

# Installing Utilities
RUN dnf install -y wget git make hostname
# Installing OpenMPI
RUN dnf install -y openmpi-devel
# Installing BLAS, LAPACK
RUN dnf install -y openblas-devel lapack-devel

# CMake installation
RUN wget -O cmake.sh https://github.com/Kitware/CMake/releases/download/v$CMAKE_VERSION/cmake-$CMAKE_VERSION-Linux-x86_64.sh
RUN sh cmake.sh --skip-license --exclude-subdir --prefix=/usr/local/

# Compilers installation
RUN dnf install -y clang gfortran

# Setting environment variables
ENV CC=/usr/bin/clang
ENV CXX=/usr/bin/clang++
ENV OMPI_CC=$CC
ENV OMPI_CXX=$CXX
ENV CC=/usr/lib64/openmpi/bin/mpicc
ENV CXX=/usr/lib64/openmpi/bin/mpic++
ENV FC=/usr/lib64/openmpi/bin/mpifort
ENV F77=/usr/lib64/openmpi/bin/mpifort
ENV F90=/usr/lib64/openmpi/bin/mpifort
ENV MPIRUNe=/usr/lib64/openmpi/bin/mpirun

# Building TPLs
RUN git clone https://github.com/Pressio/pressio-builder.git
WORKDIR /pressio-builder
RUN git checkout main
RUN sed -i 's/source .\/shared\/colors.sh/# colors commnted out/g' main_tpls.sh
RUN chmod +x main_tpls.sh
RUN ./main_tpls.sh -dryrun=no -build-mode=Release -target-dir=/usr/local -tpls=gtest,eigen,trilinos -cmake-generator-names=default,default,default

# Setting workdir to /in
WORKDIR /in