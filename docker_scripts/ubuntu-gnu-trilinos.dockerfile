ARG UBUNTU_VERSION=22.04
FROM ubuntu:${UBUNTU_VERSION}

ARG CMAKE_VERSION=3.27.8
ARG COMPILER_VERSION=11
ARG CC=gcc-${COMPILER_VERSION}
ARG CXX=g++-${COMPILER_VERSION}
ARG GFORTRAN=gfortran-${COMPILER_VERSION}

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y -q && \
    apt-get upgrade -y -q && \
    apt-get install -y -q --no-install-recommends \
        ca-certificates \
        git \
        libgtest-dev \
        liblapack-dev \
        libopenblas-dev \
        libopenmpi-dev \
        make \
        python3 \
        python3-numpy \
        wget \
        $CC $CXX $GFORTRAN && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# CMake installation
RUN wget -O cmake.sh https://github.com/Kitware/CMake/releases/download/v$CMAKE_VERSION/cmake-${CMAKE_VERSION}-Linux-x86_64.sh && \
    sh cmake.sh --skip-license --exclude-subdir --prefix=/usr/local/ && \
    rm cmake.sh

RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/${CC} 10
RUN update-alternatives --install /usr/bin/g++ g++ /usr/bin/${CXX} 10
RUN update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/${GFORTRAN} 10
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

# Building Trilinos
WORKDIR /home
RUN git clone https://github.com/trilinos/Trilinos.git && \
    cd Trilinos && \
    git checkout trilinos-release-14-4-0

RUN cmake -B Trilinos/builddir \
        -D CMAKE_BUILD_TYPE:STRING=Release \
        -D BUILD_SHARED_LIBS:BOOL=ON \
        -D TPL_FIND_SHARED_LIBS=ON \
        -D Trilinos_LINK_SEARCH_START_STATIC=OFF \
        -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
        -D TPL_ENABLE_MPI:BOOL=ON \
        -D MPI_C_COMPILER:FILEPATH=/usr/bin/mpicc \
        -D MPI_CXX_COMPILER:FILEPATH=/usr/bin/mpic++ \
        -D MPI_USE_COMPILER_WRAPPERS:BOOL=ON \
        -D Trilinos_ENABLE_Fortran:BOOL=ON \
        -D MPI_Fortran_COMPILER:FILEPATH=/usr/bin/mpifort \
        -D Trilinos_ENABLE_TESTS:BOOL=OFF \
        -D Trilinos_ENABLE_EXAMPLES:BOOL=OFF \
        -D TPL_ENABLE_BLAS=ON \
        -D TPL_ENABLE_LAPACK=ON \
        -D Kokkos_ENABLE_SERIAL:BOOL=ON \
        -D Kokkos_ENABLE_THREADS:BOOL=OFF \
        -D Kokkos_ENABLE_OPENMP:BOOL=OFF \
        -D Kokkos_ENABLE_DEPRECATED_CODE=OFF \
        -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
        -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF \
        -D Trilinos_ENABLE_Teuchos:BOOL=ON \
        -D Trilinos_ENABLE_Epetra:BOOL=ON \
        -D Trilinos_ENABLE_Tpetra:BOOL=ON \
        -D Tpetra_ENABLE_DEPRECATED_CODE:BOOL=OFF \
        -D Tpetra_ENABLE_TSQR:BOOL=ON \
        -D Trilinos_ENABLE_EpetraExt:BOOL=ON \
        -D Trilinos_ENABLE_AztecOO:BOOL=ON \
        -D Trilinos_ENABLE_Ifpack:BOOL=ON \
        -D Trilinos_ENABLE_Ifpack2:BOOL=ON \
        -D Trilinos_ENABLE_ROL:BOOL=ON \
        -D CMAKE_INSTALL_PREFIX:PATH=/home/pressio_builds/trilinos/install \
        -S Trilinos && \
    cmake --build Trilinos/builddir --target install

# Cleaning after builds
RUN rm -rf Trilinos
