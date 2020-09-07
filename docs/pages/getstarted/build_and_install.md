

# Dependencies and Installation


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

| Variable                        | Description                          | Default Value                                            |
| ------------------              | ---------------                      | -----------                                              |
| `PRESSIO_ENABLE_DEBUG_PRINT`    | to enable debugging print statements | `OFF`                                                    |
| `PRESSIO_ENABLE_EIGEN`          | self-explanatory                     | `ON`                                                     |
| `PRESSIO_ENABLE_TPL_MPI`        | self-explanatory                     | `OFF`                                                    |
| `PRESSIO_ENABLE_TPL_TRILINOS`   | self-explanatory                     | `OFF`                                                    |
| `PRESSIO_ENABLE_TPL_KOKKOS`     | self-explanatory                     | `OFF`; is set `ON` when `PRESSIO_ENABLE_TPL_TRILINOS=ON` |
| `PRESSIO_ENABLE_TEUCHOS_TIMERS` | self-explanatory                     | `OFF`                                                    |
| `PRESSIO_ENABLE_TPL_BLAS`       | self-explanatory                     | `OFF`                                                    |
| `PRESSIO_ENABLE_TPL_LAPACK`     | self-explanatory                     | `OFF`                                                    |
| `PRESSIO_ENABLE_TPL_PYBIND11`   | self-explanatory                     | `OFF`                                                    |
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


## Ok, but what are the steps to get, install and use Pressio in practice?
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
This step would process the cmake arguments and copy the `pressio`
headers to the install prefix `<where-you-cloned-pressio>/install`.
Note that nothing is built yet because `pressio` is header-only.

Since `pressio` is header-only, any TPL you want to enable
is **not** needed now during the installation of `pressio`, but it is
needed when you build any code that *uses* pressio.
<!-- * -->
<br>
<li><p>Now you point to the installed headers when building your application,
and include the `pressio` header `pressio.hpp`.
For the sake of the argument, imagine a simple main file as follows:
@m_class{m-code-figure} @parblock
@code{.cpp}
#include "pressio.hpp"
// ...
int main(){
  // do what you need
}
@endcode
</p></li>




<!-- For example, if your code is: -->
<!-- @m_class{m-code-figure} @parblock -->
<!-- @code{.cpp} -->
<!-- int main{ -->
<!--  //something that uses pressio with Trilinos enabled -->
<!-- } -->
<!-- @endcode -->
<!-- then you would first install pressio and enable the Trilinos define,  -->
<!-- clone and install Eigen and then build the code above by doing:  -->
<!-- @m_class{m-code-figure} @parblock -->
<!-- @code{.bash} -->
<!-- g++ -I  -->
<!-- @endcode -->






<!-- <li><p>Inside the `pressio/packages` subdirectory, rename the configure file as: -->
<!-- @m_class{m-code-figure} @parblock -->
<!-- @code{.bash} -->
<!-- mv pressio_cmake_config.h.in pressio_cmake_config.h -->
<!-- @endcode -->
<!-- </p></li> -->
<!-- <br> -->
<!-- <li><p>Edit the `pressio_cmake_config.h` to define the target variables you want to enable, -->
<!-- and comment all the others. -->
<!-- For example, if we want to enable all the code in `pressio` supporting Trilinos, we would do: -->
<!-- @m_class{m-code-figure} @parblock -->
<!-- @code{.cpp} -->
<!-- //#cmakedefine PRESSIO_ENABLE_DEBUG_PRINT -->
<!-- #cmakedefine PRESSIO_ENABLE_TPL_TRILINOS -->
<!-- #cmakedefine PRESSIO_ENABLE_TEUCHOS_TIMERS -->
<!-- //#cmakedefine PRESSIO_ENABLE_TPL_PYBIND11 -->
<!-- //#cmakedefine PRESSIO_ENABLE_TPL_KOKKOS -->
<!-- //#cmakedefine PRESSIO_ENABLE_TPL_BLAZE -->
<!-- //#cmakedefine PRESSIO_ENABLE_TPL_ARMADILLO -->
<!-- //#cmakedefine PRESSIO_ENABLE_TPL_BLAS -->
<!-- //#cmakedefine PRESSIO_ENABLE_TPL_LAPACK -->
<!-- //#cmakedefine PRESSIO_ENABLE_TPL_MPI -->
<!-- @endcode -->
<!-- </p></li> -->

<!-- You clone the `pressio` repo, and from your code you include -->
<!-- the `pressio/packages` to find the `pressio` headers. -->
<!-- Since `pressio` uses preprocessor directives to selectively -->
<!-- enable/disable code for target TPLs, to when you build your code you -->
<!-- need to have these preprocessor directives defined. -->

<!-- The full -->
<!-- For example, if your code uses Trilinos, to enabled the Trilinos-related code in `pressio` you need to have `PRESSIO_ENABLE_TPL_TRILINOS` defined *before* you include -->
<!-- the `pressio` headers. The list of CMake options to enable can be found [here](./list_of_cmake_optional_vars_to_enable.md). -->
