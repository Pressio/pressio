FROM intel/oneapi-hpckit:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y -q && \
    apt-get upgrade -y -q && \
    apt-get install -y -q --no-install-recommends \
        cmake \
        git \
        gnupg2 \
        libeigen3-dev \
        libgtest-dev \
        make \
        software-properties-common && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV CC=/opt/intel/oneapi/compiler/latest/linux/bin/icx
ENV CXX=/opt/intel/oneapi/compiler/latest/linux/bin/icpx

# Creating in and out directories
RUN mkdir /in
RUN mkdir /out

# Setting workdir to /in
WORKDIR /in
