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

# Building TPLs
RUN git clone https://github.com/Pressio/pressio-builder.git
WORKDIR /pressio-builder
RUN git checkout main
RUN sed -i 's/source .\/shared\/colors.sh/# colors commnted out/g' main_tpls.sh
RUN chmod +x main_tpls.sh
RUN ./main_tpls.sh -dryrun=no -build-mode=Release -target-dir=/usr/local -tpls=eigen,gtest

# Setting LD_LIBRARY_PATH variable
ENV LD_LIBRARY_PATH=/opt/intel/oneapi/vpl/2021.4.0/lib:/opt/intel/oneapi/tbb/2021.3.0/env/../lib/intel64/gcc4.8:/opt/intel/oneapi/mpi/2021.3.1//libfabric/lib:/opt/intel/oneapi/mpi/2021.3.1//lib/release:/opt/intel/oneapi/mpi/2021.3.1//lib:/opt/intel/oneapi/mkl/2021.3.0/lib/intel64:/opt/intel/oneapi/itac/2021.3.0/slib:/opt/intel/oneapi/ippcp/2021.3.0/lib/intel64:/opt/intel/oneapi/ipp/2021.3.0/lib/intel64:/opt/intel/oneapi/dnnl/2021.3.0/cpu_dpcpp_gpu_dpcpp/lib:/opt/intel/oneapi/debugger/10.1.2/gdb/intel64/lib:/opt/intel/oneapi/debugger/10.1.2/libipt/intel64/lib:/opt/intel/oneapi/debugger/10.1.2/dep/lib:/opt/intel/oneapi/dal/2021.3.0/lib/intel64:/opt/intel/oneapi/compiler/2021.3.0/linux/lib:/opt/intel/oneapi/compiler/2021.3.0/linux/lib/x64:/opt/intel/oneapi/compiler/2021.3.0/linux/lib/emu:/opt/intel/oneapi/compiler/2021.3.0/linux/lib/oclfpga/host/linux64/lib:/opt/intel/oneapi/compiler/2021.3.0/linux/lib/oclfpga/linux64/lib:/opt/intel/oneapi/compiler/2021.3.0/linux/compiler/lib/intel64_lin:/opt/intel/oneapi/ccl/2021.3.0/lib/cpu_gpu_dpcpp

# Setting workdir to /in
WORKDIR /in