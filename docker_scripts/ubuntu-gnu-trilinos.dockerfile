ARG UBUNTU_VERSION=latest
FROM ubuntu:${UBUNTU_VERSION}

ARG CC=gcc
ARG CXX=g++
ARG GFORTRAN=gfortran

RUN apt-get update -y -q && \
    apt-get upgrade -y -q && \
    apt-get install -y -q --no-install-recommends \
        ca-certificates \
        cmake \
        git \
        libeigen3-dev \
        libgtest-dev \
        liblapack-dev \
        libopenblas-dev \
        make \
        openmpi-bin \
        openmpi-doc \
        python3 \
        python3-numpy \
        wget \
        $CC $CXX $GFORTRAN && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

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
RUN ./main_tpls.sh -dryrun=no -build-mode=Release -target-dir=../../pressio_builds -tpls=trilinos -cmake-generator-names=default
RUN git reset --hard origin/main

# Cleaning after builds
WORKDIR /home
# RUN rm -rf pressio_builds/trilinos/build && rm -rf pressio_builds/trilinos/Trilinos

# Setting workdir to /
WORKDIR /