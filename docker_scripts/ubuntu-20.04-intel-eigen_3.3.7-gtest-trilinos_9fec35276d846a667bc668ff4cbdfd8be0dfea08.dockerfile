FROM ubuntu:focal
# Declaring build variables
ARG TZ=Europe/Warsaw
ARG CMAKE_VERSION=3.16.8
ARG DEBIAN_FRONTEND=noninteractive

# Setting timezone
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Creating in and out directories
RUN mkdir /in
RUN mkdir /out

# System update and packages installation
RUN apt-get update && apt-get upgrade -y
# Installing Utilities
RUN apt-get install -y wget git make gnupg2 software-properties-common apt-utils
# Installing OpenMPI
RUN apt-get install -y openmpi-bin openmpi-doc

# CMake installation
RUN wget -O cmake.sh https://github.com/Kitware/CMake/releases/download/v$CMAKE_VERSION/cmake-$CMAKE_VERSION-Linux-x86_64.sh
RUN sh cmake.sh --skip-license --exclude-subdir --prefix=/usr/local/

# Compilers installation
RUN wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN add-apt-repository "deb https://apt.repos.intel.com/oneapi all main"
RUN apt-get install -y intel-hpckit

# Intel MKL installation and configuration
RUN add-apt-repository "deb https://apt.repos.intel.com/mkl all main"
RUN apt-get install -y -qq intel-mkl-2020.4-912

# Intel Blas and Lapack configuration
RUN update-alternatives --install /usr/lib/x86_64-linux-gnu/libblas.so libblas.so-x86_64-linux-gnu /opt/intel/mkl/lib/intel64/libmkl_rt.so 50
RUN update-alternatives --install /usr/lib/x86_64-linux-gnu/libblas.so.3 libblas.so.3-x86_64-linux-gnu /opt/intel/mkl/lib/intel64/libmkl_rt.so 50
RUN update-alternatives --install /usr/lib/x86_64-linux-gnu/liblapack.so liblapack.so-x86_64-linux-gnu /opt/intel/mkl/lib/intel64/libmkl_rt.so 50
RUN update-alternatives --install /usr/lib/x86_64-linux-gnu/liblapack.so.3 liblapack.so.3-x86_64-linux-gnu /opt/intel/mkl/lib/intel64/libmkl_rt.so 50
RUN echo "/opt/intel/lib/intel64" > /etc/ld.so.conf.d/mkl.conf
RUN echo "/opt/intel/mkl/lib/intel64" >> /etc/ld.so.conf.d/mkl.conf
RUN ldconfig

# Setting environment variables
ENV MKLROOT=/opt/intel/compilers_and_libraries_2020.4.304/linux/mkl
ENV LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2020.4.304/linux/tbb/lib/intel64_lin/gcc4.7:/opt/intel/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64_lin:/opt/intel/compilers_and_libraries_2020.4.304/linux/mkl/lib/intel64_lin
ENV LIBRARY_PATH=/opt/intel/compilers_and_libraries_2020.4.304/linux/tbb/lib/intel64_lin/gcc4.7:/opt/intel/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64_lin:/opt/intel/compilers_and_libraries_2020.4.304/linux/mkl/lib/intel64_lin
ENV NLSPATH=/opt/intel/compilers_and_libraries_2020.4.304/linux/mkl/lib/intel64_lin/locale/%l_%t/%N
ENV CPATH=/opt/intel/compilers_and_libraries_2020.4.304/linux/mkl/include
ENV PKG_CONFIG_PATH=/opt/intel/compilers_and_libraries_2020.4.304/linux/mkl/bin/pkgconfig
ENV CC=/opt/intel/oneapi/compiler/latest/linux/bin/intel64/icc
ENV CXX=/opt/intel/oneapi/compiler/latest/linux/bin/intel64/icpc
ENV OMPI_CC=$CC
ENV OMPI_CXX=$CXX
ENV CC=/usr/bin/mpicc
ENV CXX=/usr/bin/mpic++
ENV FC=/usr/bin/mpifort
ENV F77=/usr/bin/mpifort
ENV F90=/usr/bin/mpifort
ENV MPIRUNe=/usr/bin/mpirun

# Building TPLs
RUN git clone https://github.com/Pressio/pressio-builder.git
WORKDIR /pressio-builder
RUN git checkout main
RUN sed -i 's/source .\/shared\/colors.sh/# colors commnted out/g' main_tpls.sh
RUN chmod +x main_tpls.sh
RUN ./main_tpls.sh -dryrun=no -build-mode=Release -target-dir=/usr/local -tpls=gtest,eigen,trilinos -cmake-generator-names=default,default,default

# Setting workdir to /in
WORKDIR /in