

# Installation and Dependencies


@m_class{m-block m-info}

@parblock This page describes the dependencies of `pressio` and its installation process.
By the end, you should be able to get pressio, install it,
and point to the installed headers from your application.
@endparblock


@m_class{m-block m-info}

@parblock `pressio` is header-only, so to use it one does not need to precompile it
into a library and linking to it.
However, since we use preprocessor directives to conditionally
enable/disable code based on target third-party libraries,
one needs to account for this. See below for the details.
@endparblock



## Dependencies

Some packages of `pressio` contain code and implementations
that are specific to third-party libraries (TPLs).
For example, the `containers` and `ops` packages of `pressio` contain
thin wrappers and kernels that are custom-built for target libraries.
The main reason for doing this is that we aim to alleviate the user from
writing custom operations and allows `pressio` to decide when and how to leverage
the native libraries' operations to obtain the best performance.
This should facilitate the integration and use of `pressio` by existing applications.
Obviously, this is a growing capability and we currently only
provide built-in support to some external HPC libraries (see below).
We can distinguish between *optional* and *required* dependencies.

| TPL Library Name   | Optional/Required | Version Known to Work |
| ------------------ | ---------------   |                       |
| Eigen              | Required          | 3.3.7                 |
| Trilinos           | Optional          | 12.17.00              |
| MPI                | Optional          | --                    |
| Kokkos             | Optional          | 3.1.0                 |
| BLAS               | Optional          | --                    |
| LAPACK             | Optional          | --                    |
| Pybind11           | Optional          | v2.5                  |
| GoogleTest         | Optional          | 1.10.0                |
|                    |                   |                       |


Enabling/disabling specific dependencies is done via the following cmake variables:

| Variable                        | Description                          | Default Value                                                                                   |
| ------------------              | ---------------                      | -----------                                                                                     |
| `PRESSIO_ENABLE_TPL_EIGEN`      | self-explanatory                     | `ON`                                                                                            |
| `PRESSIO_ENABLE_TPL_TRILINOS`   | self-explanatory                     | `OFF`                                                                                           |
| `PRESSIO_ENABLE_TPL_MPI`        | self-explanatory                     | `OFF`  automatically `ON` if `PRESSIO_ENABLE_TPL_TRILINOS=ON`                                   |
| `PRESSIO_ENABLE_TPL_KOKKOS`     | self-explanatory                     | `OFF`; automatically `ON` if `PRESSIO_ENABLE_TPL_TRILINOS=ON`                                   |
| `PRESSIO_ENABLE_TEUCHOS_TIMERS` | self-explanatory                     | `OFF`  automatically `ON` if `PRESSIO_ENABLE_TPL_TRILINOS=ON`                                   |
| `PRESSIO_ENABLE_TPL_BLAS`       | self-explanatory                     | `OFF`; automatically `ON` if `PRESSIO_ENABLE_TPL_LAPACK=ON` or `PRESSIO_ENABLE_TPL_TRILINOS=ON` |
| `PRESSIO_ENABLE_TPL_LAPACK`     | self-explanatory                     | `OFF`; automatically `ON` if `PRESSIO_ENABLE_TPL_BLAS=ON` or `PRESSIO_ENABLE_TPL_TRILINOS=ON`   |
| `PRESSIO_ENABLE_TPL_PYBIND11`   | self-explanatory                     | `OFF`                                                                                           |
| `PRESSIO_ENABLE_DEBUG_PRINT`    | to enable debugging print statements | `OFF`                                                                                           |
<!-- | `PRESSIO_ENABLE_TPL_TORCH`| self-explanatory | `OFF` |-->

@m_class{m-block m-default}

@par
	Eigen is the only required dependency because it is the
	default choice for instantiating the ROM data structures
	and solving the (dense) ROM problem.



Obviously, the choice of which TPLs to enable is related to
your application's dependency requirements.
For example, if you have an application that relies on
Trilinos data structures and want to use `pressio`,
then it makes sense to enable the Trilinos dependency.
On the contrary, if you have an application that relies only on
Eigen data structures, then it makes sense to only leave only Eigen on
and disable the rest.

Also, we note that some of the cmake variables listed above are connected
and cannot be turned on individualy.
For example, if we enable Trilinos then `pressio` automatically
enables also Kokkos, BLAS, LAPACK and MPI.
The reason for this choice is that in a production scenario---which is what
pressio mainly targets---it is reasonable
to expect a user to have Trilinos built with BLAS, LAPACK, MPI and Kokkos support.

There might be other constraints on the variables one can set.
The reason for this is twofold: (a) to simplify what the user needs
to provide; and (b) we belive some of these constraints are necessary, like
the Trilinos example above or always requiring BLAS and LAPACK to be simulateneously on.
<!-- Note that, since `pressio` is header-only, any TPL you want to enable -->
<!-- is not really needed when installing `pressio`, but it is needed when -->
<!-- you build any code that *uses* pressio. -->
<!-- Therefore, you need to make sure that before you use `pressio` in your code, -->
<!-- you include/link to any TPL you want to use. -->
<!-- At the very minimum, you need to have Eigen installed. -->


## In practice, what are the steps to get, install and use Pressio?
We suggest to follow these steps:
<ol>
<li>Clone [pressio](https://github.com/Pressio/pressio) (defaults to the master branch)</li>

<li>Create a build and install subdirectories
@m_class{m-code-figure}

@code{.bash}
cd <where-you-cloned-pressio>
mkdir build && mkdir install
@endcode
</li>

<li> Use cmake to configure by passing to the comand line the target
list of cmake variables to define. For example, if we want to enable
support in `pressio` for Trilinos and the debug prints, we would do:
@m_class{m-code-figure}

@code{.bash}
export PRESSIO_SRC=<where-you-cloned-pressio>
cd <where-you-cloned-pressio>/build

cmake -D CMAKE_INSTALL_PREFIX=../install \
	  -D PRESSIO_ENABLE_TPL_TRILINOS=ON \
	  -D PRESSIO_ENABLE_DEBUG_PRINT=ON \
	  ${PRESSIO_SRC}

make install # nothing is built, just headers copied to installation
@endcode
</li>

Note that this step does **not** build anything because `pressio` is header-only,
but only processes the cmake arguments and copies the pressio headers to the
install prefix `<where-you-cloned-pressio>/install`.
If you want, you can inspect the file `<where-you-cloned-pressio>/install/presssio_cmake_config.h`
which contains the cmake variables defined.

We also remark that during the step above pressio does not need to know
if and where a target TPL exists in your system.
Above you are simply telling Pressio that you have
a certain TPL and want to enable the corresponding code in pressio for later use.
Those TPLs will be needed when you build any code that *uses* pressio.

<li> When building your application, you point to the installed headers
and include the `pressio` header `pressio.hpp`, for example as:
@m_class{m-code-figure}

@code{.cpp}
#include "pressio.hpp"
// ...
int main(){
  // do what you need
}
@endcode
</li>
</ol>


@m_class{m-block m-warning}

@par Warning:
The procedue above is highly advised because it enables `pressio`
to properly process the cmake options and turn on/off based
on certain conditions (as explained above).
The alternative way to use pressio would be to just clone the repo,
point to its source code and use cmake variables directly when building
your code. However, this could have unexpected consequences since
you would be resposible to set the variables correctly but you would not
know exactly all the possible constrants.
Therefore, we suggest to use the steps above.
