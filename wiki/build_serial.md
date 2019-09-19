
# Building and Running pressio-tutorials
This step-by-step page explains how to build the tutorials from scratch on Unix.

*Disclaimer*: the guide below does **not** assume you are a Unix/CS/coding ninja, rather the opposite. It is written with the goal to provide a complete and detailed guide without taking any step for granted. As such, if you are a Unix/CS/coding ninja, some steps will be obvious to you, so please do not get offended!

## Prerequisites
In order for the steps below to be successful, you need:

* C/C++ compilers: either Clang or GNU.
The current version of the tutorials does NOT need MPI. But if you have MPI compiler wrappers, you can use those to build.

* CMake, minimun version 2.8.12.

**Note**: we purposefully avoid (for the time being) using distributed data structures in these tutorials because we want a simpler framework to explain some concepts. However, we remark that the code shown in the tutorials can be used for MPI-based distributed data structures with almost no change. This will be explained in more detail within some of the tutorials.


<!---------------------------------------------------->
## 1. Prep
<!-- For the sake of clarity, let us assume your name is *John Doe*, and your username is `johndoe`. -->

(a) Create (or choose) a directory where you want to clone all repos needed and where to put all builds, e,g.:
```bash
mkdir $HOME/pressio_repos
mkdir $HOME/pressio_builds
```

(b) To make things easier and cleaner below, create environment variables to refer to these directories:
```bash
export PRESSIO_REPOS=$HOME/pressio_repos
export PRESSIO_BUILDS=$HOME/pressio_builds
```

(c) Unless you already have them, set the following compilers environment variable:
```bash
export CC=<path-to-your-C-compiler>
export CXX=<path-to-your-CXX-compiler>
```
These are needed because `CC` and `CXX` are used to do all the builds.


<!---------------------------------------------------->
## 2. Cloning

To build the tutorials, you need to clone the following repos:
```bash
cd ${PRESSIO_REPOS}
git clone git@github.com:Pressio/pressio-builder.git
git clone --recursive git@github.com:Pressio/pressio.git
git clone git@github.com:Pressio/pressio-tutorials.git
```

<!---------------------------------------------------->
## 3. Install TPLs

We only need Eigen (for now), so you can simply run the command:
```bash
cd ${PRESSIO_REPOS}/pressio-builder
./main_tpls.sh -dryrun=0 -tpls=eigen -target-dir=${PRESSIO_BUILDS}
```
To learn more about the script's command line args, type `./main_tpls.sh -h`.

<!---------------------------------------------------->
## 4. Install Pressio
From the same directory, i.e. `${PRESSIO_REPOS}/pressio-builder`, run the command:
```bash
./main_pressio.sh \
	-dryrun=0 \
	-pressio-src=${PRESSIO_REPOS}/pressio \
	-packages=rom \
	-target-dir=${PRESSIO_BUILDS} \
	-eigen-path=${PRESSIO_BUILDS}/eigen/install \
	-cmake-generator-name=default_for_tutorials
```
To learn more about the script's command line args, type `./main_pressio.sh -h`.

<!---------------------------------------------------->
## 5. Building the tutorials
From the same directory, i.e. `${PRESSIO_REPOS}/pressio-builder`, run the command:
```bash
./main_tutorials.sh \
	-dryrun=0 \
	-pressio-tutorials-src=${PRESSIO_REPOS}/pressio-tutorials \
	-target-dir=$HOME/Desktop/pressio_builds \
	-build-mode=Release \
	-eigen-path=${PRESSIO_BUILDS}/eigen/install \
	-pressio-path=${PRESSIO_BUILDS}/pressio/install
```
To learn more about the script's command line args, type `./main_tutorials.sh -h`.
