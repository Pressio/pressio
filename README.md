
# Overview

Pressio is a collection of repositories providing reduced-order models (ROMs) capabilties.
Specifically, the whole Pressio framework currently includes the following repositories:

* `pressio`: the main C++ repo containing the ROM algorithms and supporting functionalities;

* `pressio-tutorials`: C++ tutorials explaining how to use `pressio` and its functionalities;

* `pressio-builder`: an auxiliary repo with bash helper scripts for configuring/building/installing `pressio`, and `pressio-tutorials`.

======================================================================================
## License
Pressio is released with the following [LICENSE](./LICENSE).

======================================================================================
## Cloning
To clone `pressio` (provided you have Git) use: `git clone https://github.com/Pressio/pressio.git`

======================================================================================
## Structure
For a description of `pressio` code structure, see [here](https://github.com/Pressio/pressio/wiki/Structure-of-pressio).

======================================================================================
## Installing
`pressio` is a header-only library, so there is **no building process needed** if you need to use it from your code. 
You just clone the code, and point to the headers which are inside the packages directory. 
However, since `pressio` uses preprocessor directives to selectively enable/disable code for target TPLs, when you build your code you need to account for this. Also, if you want to build the tests in `pressio`, in that case too you need to enable specific TPLs if you want related tests to be turned on. 

For a list of CMake options to enable see [this file](./list_of_cmake_optional_vars_to_enable.md).

Sample cmake configure lines can be found [here](https://github.com/Pressio/pressio/wiki/Sample-CMake-configure-lines-for-pressio).
<!--
======================================================================================
## TPLs
At the time of this writing, `pressio` has only one required dependency, namely Eigen, and a few **optional** ones, namely Gtest, Pybind11, Trilinos, Kokkos, Armadillo, Blas, Lapack, Blaze. This set of TPLs will liekly increase over time as we add support for more external packages, e.g. Petsc. However, one of the key choices is that we will keep most dependencies optional. Moreover, these TPLs are NOT needed for the installation process since `pressio` is header-only.

======================================================================================
## Configuring and Building
Configuring and building `pressio` can be done in two ways:

* since `pressio` uses CMake, you can use a typical CMake configure/build/install process. Note that some TPLs are needed.

* by clonig and using the following helper repo: `git clone https://github.com/Pressio/pressio-builder.git`
The advantage of using the helper repo is that it automates the installation of TPLs.

### Sample building steps

This wiki will be updated over time, but to get started, we provide here basic references to configure `pressio` and build its tests for a few target scenarios.

* Follow [this](./wiki/build_serial.md) for a basic *serial* build that uses only Gtest and Eigen and it is done with `pressio-builder` (which automatically builds) Gtest, Eigen for you.
 -->

======================================================================================
## Disclaimer

* Pressio is work-in-progress. At the time of this writing, it is a fairly young project and things are obviously evolving. Several package would benefit from substantial work on testing and documentation, and this is ongoing. However, `pressio` is functional and has been already tested/deployed on large-scale applications.

*This document is in progress, more details soon.*
