
# Build Tests with Eigen only


@m_class{m-block m-info}

@par What does this page describe?
This page shows the steps to build the unit and regression tests
in `pressio` depending only on Eigen (and GoogleTest).
By the end, you should be able to have pressio built
with GoogleTest and Eigen enabled, run the test suite inside.



@m_class{m-block m-info}

@par Why can this be useful?
If you just want a quick way to get the code up-and-running and
ready to play around with, these steps are what you want.
Once you have the build ready, you can easily test/play/explore.



@m_class{m-block m-warning}

@par Disclaimer
The guide below does **not** assume you are
a Unix/CS/coding ninja, rather the opposite. It is written with the goal
to provide a detailed guide without taking any step for granted.
As such, if you are a Unix/CS/coding ninja, some steps will be fairly obvious to you!



## Prerequisites

In order for the steps below to be successful, you need:

* C and C++11 compilers: either Clang or GNU
* CMake >= 3.11.0
* Bash >= 3.2.57


## 1. Prep

Create (or choose) a directory where you want to clone all repos needed and where to put all builds, e,g.:
```bash
mkdir $HOME/pressio_repos
mkdir $HOME/pressio_builds
```

To make things easier and cleaner below, create environment variables to refer to these directories:
```bash
export PRESSIO_REPOS=$HOME/pressio_repos
export PRESSIO_BUILDS=$HOME/pressio_builds
```

Unless you already have them, set the following compilers environment variable:
```bash
export CC=<path-to-your-C-compiler>
export CXX=<path-to-your-CXX-compiler>
```
These are needed because `CC` and `CXX` are used to do all the builds.


## 2. Cloning

You need to clone the following repos:
```bash
cd ${PRESSIO_REPOS}
git clone git@github.com:Pressio/pressio-builder.git
git clone git@github.com:Pressio/pressio.git
```

## 3. Install TPLs

To get and installed Eigen and Gtest, you can simply run the command:
```bash
cd ${PRESSIO_REPOS}/pressio-builder
./main_tpls.sh -dryrun=no -tpls=eigen,gtest -target-dir=${PRESSIO_BUILDS}
```
<!-- To learn more about the script's command line args, type `./main_tpls.sh -h`. -->


## 4. Build the tests

From the same directory, i.e. `${PRESSIO_REPOS}/pressio-builder`, run the command:
```bash
./main_pressio.sh \
	-dryrun=no \
	-pressio-src=${PRESSIO_REPOS}/pressio \
	-target-dir=${PRESSIO_BUILDS} \
	-gtest-path=${PRESSIO_BUILDS}/gtest/install \
	-eigen-path=${PRESSIO_BUILDS}/eigen/install
	-cmake-generator-name=default_with_tests
```
<!-- To learn more about the script's command line args, type `./main_pressio.sh -h`. -->


## 5. Running the tests
To run the tests, you can follow this:
```bash
cd ${PRESSIO_BUILDS}/pressio/build
ctest
```
To learn more about ctest, you can do `ctest --help`.
