
# Pressio

Pressio is a collection of repositories providing reduced-order models (ROMs) capabilties.
Specifically, the whole Pressio framework currently includes the following *public* repositories: 

* `pressio`: this is the main C++ repo containing all the algorithms;

* `pressio-tutorials`: C++ tutorials explaining how to use `pressio` and its functionalities; 

* `pressio-builder`: an auxiliary repo with bash scripts to help with configuration/building/installation of `pressio`, and `pressio-tutorials`.


## License
Pressio is released with the following [LICENSE](./LICENSE).

## Cloning
To clone `pressio`, use this command:
```
git clone --recursive https://github.com/Pressio/pressio.git
```
The recursive option is necessary to clone a git submodule for TriBITS.
TriBITS (https://tribits.org/) provides the development framework for Pressio.

## Structure 

`pressio` is the main code repository and currently includes the following packages: 

* `mpl`: metaprogramming functionalities;

* `utils`: common functionalities, e.g., I/O helpers, static constants, etc;

* `containers`: data strutures wrappers and linear algebra;

* `apps`: suites of mini-apps used for basic testing;

* `qr`: QR factorization functionalities;

* `svd`: singular value decomposition (SVD) functionalities;

* `solvers`: linear and non-linear solvers (e.g. Gauss-Newton with and without line-search, etc.);

* `ode`: explicit and implict time steppers and integrators;

* `rom`: reduced-order modeling algorithms.

The top-down order used above is informative of the packages' dependency structure and mutual visibility. For example, every package depends on `mpl`,
but `qr` depends only on `mpl`, `utils`, `containers`. At the bottom of the hierarchy we have the `rom` package which requires all the others.
Each package contains corresponding unit- and regular tests. Splitting the framework into separate packages has several benefits.
* Maintability: `pressio` can be more easily developed and maintained since packages depend on one another through well-defined public interfaces, and appropriate namespaces are used to organize classes within each package.

* Selective usability: This modular framework allows users, if needed, to leverage invidual packages (similarly to Trilinos). For instance, if a user needs/wants just the QR methods, they simply use that package, and all the dependencies are enabled automatically.



## Configuring and Building
Configuring and building Pressio can be done in two ways: 

* by clonig and using the following helper repo: `git clone https://github.com/Pressio/pressio-builder.git`

* since Pressio uses TriBITS (which uses CMake), you can use a typical CMake configure/build/install process. Note that some TPLs are needed. 

The advantage of using the helper repo is that it automates the installation of TPLs. 

### Sample building steps

This wiki will be updated over time, but to get started, we provide here basic references 
to configure `pressio` and build its tests for target scenarios. 

* For a serial build using only Gtest, Eigen: [here](./wiki/serial_build.md)



## Disclaimer

* Pressio is work-in-progress. At the time of this writing, it is a fairly young project and things are obviously evolving. Several package would benefit from substantial work on testing and documentation, and this is ongoing. However, `pressio` is functional and has been already tested/deployed on large-scale applications.

*This document is in progress, more detailed info to come soon.*







