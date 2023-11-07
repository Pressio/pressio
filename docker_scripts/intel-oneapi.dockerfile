FROM intel/oneapi-hpckit:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | apt-key add -

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
