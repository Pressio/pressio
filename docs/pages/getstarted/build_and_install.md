

# Installation and Dependencies


@m_class{m-block m-info}

@par What does this page cover?
This page describes the dependencies of `pressio` and its installation process.
By the end, you should be able to get pressio, install it,
and point to the installed headers from your application.
@endparblock


@m_class{m-block m-info}

@par Overview
`pressio` is header-only, so to use it one does not need to precompile it
into a library and linking to it.
However, since we use preprocessor directives to conditionally
enable/disable code based on target third-party libraries,
one needs to account for this. See below for the details.
@endparblock


## Dependencies

`pressio` contains code that is specific to third-party libraries (TPLs)
to benefit the integration into existing applications.
We can distinguish between *optional* and *required* dependencies.
<!-- This is important, because it alleviates the user from writing custom operations -->
<!-- and allows `pressio` to decide when to leverage the native libraries' operations to -->
<!-- obtain the best performance. -->

| TPL Library Name   | Optional/Required |
| ------------------ | ---------------   |
| Eigen              | Required          |
| Trilinos           | Optional          |
| MPI                | Optional          |
| Kokkos             | Optional          |
| BLAS               | Optional          |
| LAPACK             | Optional          |
| Pybind11           | Optional          |

@m_class{m-block m-default}

@par
	Eigen is the only required dependency because it is the
	default choice for instantiating the (dense) ROM operators
	and solving the (dense) ROM problem.
@endparblock

Enabling/disabling specific dependencies is done via the following cmake variables:
\todo fix the autoamtic setting of pressio based on what is selected.

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

Obviously, the choice of TPLs to enable is related to the
dependency requirements of your application.
For example, if you have an application that relies on
Trilinos data structures and want to use `pressio`, then it makes sense
to enable the Trilinos dependency.
If you have an application that relies on Eigen data structures,
then it makes sense to only leave only Eigen on.

Also, some of the cmake variables listed above are connected.
For example, if we enable Trilinos then `pressio` automatically
enables also Kokkos, BLAS, LAPACK and MPI.
The reason for this is twofold: (a) to simplify what the user needs
to provide; and (b) we belive some of these constraints are necessary.



<!-- Note that, since `pressio` is header-only, any TPL you want to enable -->
<!-- is not really needed when installing `pressio`, but it is needed when -->
<!-- you build any code that *uses* pressio. -->
<!-- Therefore, you need to make sure that before you use `pressio` in your code, -->
<!-- you include/link to any TPL you want to use. -->
<!-- At the very minimum, you need to have Eigen installed. -->


## In practice, what are the steps to get, install and use Pressio?
We suggest to use the following steps to install the `pressio` headers
such that you can then point to it from within your code:
<ol>
<li><p>Clone `pressio`</p></li>
<!-- * -->
<li><p>Create a build and install subdirectories</p></li>
@m_class{m-code-figure} @parblock
@code{.bash}
cd <where-you-cloned-pressio>
mkdir build && mkdir install
@endcode
</p></li>
<!-- * -->
<br>
<li><p>Use cmake to configure by passing to the comand line the target
list of cmake variables to define. For example, if we want to enable
support in `pressio` for Trilinos, one would do:</p></li>
@m_class{m-code-figure} @parblock
@code{.bash}
cd <where-you-cloned-pressio>/build
cmake -D CMAKE_INSTALL_PREFIX=../install \
	  -D PRESSIO_ENABLE_TPL_TRILINOS=ON ..
make install # nothing is built, just headers copied to installation
@endcode
</p></li>
Note that during this step **nothing** is built because `pressio` is header-only.
This step only processes the cmake arguments and copy the headers to the
install prefix `<where-you-cloned-pressio>/install`, where you can
also find a `presssio_cmake_config.h` file containing the cmake variables defined.

Since `pressio` is header-only, any TPL you want to enable
is **not** really needed during this installation of `pressio`,
but obviously will be needed when you build any code that *uses* pressio and that TPL.
<!-- * -->
<br>
<li><p>Now you point to the installed headers when building your application,
and include the `pressio` header `pressio.hpp`, for example as:
@m_class{m-code-figure} @parblock
@code{.cpp}
#include "pressio.hpp"
// ...
int main(){
  // do what you need
}
@endcode
</p></li>
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
@endparblock
