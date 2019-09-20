
# Building andr unning serially some target tests in `pressio`
This step-by-step page explains how to build the ROM tests in `pressio` from scratch on Unix.

*Disclaimer*: the guide below does **not** assume you are a Unix/CS/coding ninja, rather the opposite. It is written with the goal to provide a complete and detailed guide without taking any step for granted. As such, if you are a Unix/CS/coding ninja, some steps will be obvious to you, so please do not get offended!

## Prerequisites
In order for the steps below to be successful, you need:

* C/C++ compilers: either Clang or GNU.
The current version of the tutorials does NOT need MPI. But if you have MPI compiler wrappers, you can use those to build.

* CMake, minimun version 3.11.0


## 1. Prep

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

You need to clone the following repos:
```bash
cd ${PRESSIO_REPOS}
git clone git@github.com:Pressio/pressio-builder.git
git clone --recursive git@github.com:Pressio/pressio.git
```

<!---------------------------------------------------->
## 3. Install TPLs

We only need Eigen and Gtest (for now), so you can simply run the command:
```bash
cd ${PRESSIO_REPOS}/pressio-builder
./main_tpls.sh -dryrun=no -tpls=eigen,gtest -target-dir=${PRESSIO_BUILDS}
```
To learn more about the script's command line args, type `./main_tpls.sh -h`.

<!---------------------------------------------------->
## 4. Build Tests in `pressio`
From the same directory, i.e. `${PRESSIO_REPOS}/pressio-builder`, run the command:
```bash
./main_pressio.sh \
	-dryrun=no \
	-pressio-src=${PRESSIO_REPOS}/pressio \
	-target-dir=${PRESSIO_BUILDS} \
	-gtest-path=${PRESSIO_BUILDS}/gtest/install \
	-eigen-path=${PRESSIO_BUILDS}/eigen/install \
	-cmake-generator-name=default \
	-package-name=rom
```
To learn more about the script's command line args, type `./main_pressio.sh -h`.

<!---------------------------------------------------->
## 5. Running the tests in `pressio`
To run the tests, you can follow this: 
```bash
${PRESSIO_BUILDS}/pressio/build
ctest 
```



