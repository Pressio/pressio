
# Installation and Dependencies


@m_class{m-block m-success}

@parblock `pressio` is header-only, so it does not need to be precompiled and linked to.
However, since we use preprocessor directives to conditionally
enable/disable code based on target third-party libraries,
one needs to account for this. See below for the details.
@endparblock



## Dependencies

Some parts of `pressio` contain code and implementations
that are specific to third-party libraries (TPLs).
An example is `pressio/ops`, which contains kernels specializations
for widely-used HPC libraries (e.g. Trilinos, Kokkos).
The main reason for doing this is that we aim, where possible,
to alleviate the user from writing custom operations and allow `pressio` to decide when and how to leverage
the native libraries' operations to obtain the best performance.
This should facilitate the integration and use of `pressio` by existing applications.
This is a growing capability and we currently only
provide built-in support to some external HPC libraries (see below).
Obviously, these TPLs-specific specializations need to be guarded with
preprecessor directives, and enabled only if one can access the TPLs.

Enabling/disabling specific dependencies can be done via the following cmake variables:

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

### Optional vs Required

Eigen is the only required dependency because it is the
default choice for instantiating the ROM data structures
and solving the (dense) ROM problem.

| TPL Library Name   | Optional/Required | Version Known to Work |
| ------------------ | ---------------   |                       |
| Eigen              | Required          | 3.3.7                 |
| Trilinos           | Optional          | 12.17.00              |
| MPI                | Optional          | --                    |
| Kokkos             | Optional          | 3.1.0                 |
| BLAS               | Optional          | --                    |
| LAPACK             | Optional          | --                    |
| Pybind11           | Optional          | v2.6                  |
| GoogleTest         | Optional          | 1.10.0                |
|                    |                   |                       |


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
<ol>
<li>Clone [pressio](https://github.com/Pressio/pressio) (defaults to the main branch),
or you can pick a [release version](https://github.com/Pressio/pressio/releases) </li>

<li>Create a build and install subdirectories
@m_class{m-code-figure}

@code{.bash}
cd <where-you-cloned-pressio>
mkdir build && mkdir install
@endcode
</li>

<li> Use cmake to configure by passing to the comand line the target list of cmake variables to define. <br/>
For example, suppose we want to enable support for Trilinos and the debug prints. We would do:
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

*Note that this step does **not** build anything because `pressio` is header-only,
but only processes the cmake arguments and copies the pressio headers to the
install prefix* `<where-you-cloned-pressio>/install`.<br/>
If you want, inspect the file `<where-you-cloned-pressio>/install/presssio_cmake_config.h`
which contains the cmake variables configuration.

We also remark that during the configuration step above pressio
does not need to know where a target TPL exists in your system.
In the configuration step above, you are simply telling Pressio that you have
a certain TPL and want to enable the corresponding code in pressio.
The TPLs will be needed at linking stage when you build *your* code that *uses* pressio.

<li> When building your application to use pressio, you just have to point to
the installation directory `<where-you-cloned-pressio>/install` with the installed
pressio headers, and you can access all pressio functionalities via the C++ include `#include<pressio.hpp>`:
@m_class{m-code-figure}

@code{.cpp}
#include "pressio/what_you_need.hpp"
// ...
int main(){
  // do something
}
@endcode
</li>
</ol>

@m_class{m-block m-warning}

@par
The procedure above is advised because it enables `pressio`
to properly process the cmake options and turn on/off based
on certain conditions (as explained above).
The alternative way to use pressio would be to just clone the repo,
point to the `<where-you-cloned-pressio>/include` subdirectory
and use cmake variables directly when building your code.
However, this could have unexpected consequences since
you would be resposible to set the variables correctly but you would not
know exactly all the possible constraints.
