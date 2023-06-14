.. role:: raw-html-m2r(raw)
   :format: html

Installation and Dependencies
=============================

.. tip::

    pressio is header-only, so it does not need to be precompiled and linked to.
    However, since we use preprocessor directives to conditionally
    enable/disable code based on target third-party libraries,
    one needs to account for this. See below for the details.

.. warning::

    To use pressio, you need at least a C++17 compiler.

Dependencies
------------

Some parts of ``pressio`` contain code and implementations
that are specific to third-party libraries (TPLs).
An example is ``pressio/ops``\ , which contains kernels specializations
for widely-used HPC libraries (e.g. Trilinos, Kokkos).
The main reason for doing this is that we aim, where possible,
to alleviate the user from writing custom operations and allow ``pressio`` to decide when and how to leverage
the native libraries' operations to obtain the best performance.
This should facilitate the integration and use of ``pressio`` by existing applications.
This is a growing capability and we currently only
provide built-in support to some external HPC libraries (see below).
Obviously, these TPLs-specific specializations need to be guarded with
preprecessor directives, and enabled only if one can access the TPLs.

Enabling/disabling specific dependencies can be done via
the cmake variables `listed here <keywords.html>`__.


Optional vs Required
^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 10 50 40
   :align: left

   * - TPL Name
     - Optional/Required
     - Version Known to Work/run in CI
   * - Eigen
     - Required
     - 3.3.7
   * - Trilinos
     - Optional
     - commit: ef73d14babf6e7556b0420add98cce257ccaa56b
   * - MPI
     - Optional
     - --
   * - Kokkos
     - Optional
     - 3.1.0
   * - BLAS
     - Optional
     - --
   * - LAPACK
     - Optional
     - --
   * - GoogleTest
     - Optional
     - 1.10.0

Eigen is the only required dependency because it is the
default choice for instantiating the ROM data structures
and solving the (dense) ROM problem.

In practice, what are the steps to get, install and use Pressio?
----------------------------------------------------------------

\todo: add description for using pressio without configuring,
so one has to define the options directly when configuring
their code or inside the code directly.

1. Clone `pressio <https://github.com/Pressio/pressio>`_ (defaults to the main branch),
or you can pick a `release version <https://github.com/Pressio/pressio/releases>`_

2. Create a build and install subdirectories

.. code-block:: bash

   cd <where-you-cloned-pressio>
   mkdir build && mkdir install

3. Use cmake to configure by passing to the comand line the target list of cmake variables to define.

For example, suppose we want to enable support for Trilinos and the debug prints. We would do:

.. code-block:: bash

   export PRESSIO_SRC=<where-you-cloned-pressio>
   cd <where-you-cloned-pressio>/build

   cmake -D CMAKE_INSTALL_PREFIX=../install \
         -D PRESSIO_ENABLE_TPL_TRILINOS=ON \
         -D PRESSIO_ENABLE_DEBUG_PRINT=ON \
         ${PRESSIO_SRC}

   make install # nothing is built, just headers copied to installation

Note that this step does **not** build anything because ``pressio`` is header-only,
but only processes the cmake arguments and copies the pressio headers to the
install prefix ``<where-you-cloned-pressio>/install``.
If you want, inspect the file ``<where-you-cloned-pressio>/install/presssio_cmake_config.h``
which contains the cmake variables configuration.

We also remark that during the configuration step above pressio
does not need to know where a target TPL exists in your system.
In the configuration step above, you are simply telling Pressio that you have
a certain TPL and want to enable the corresponding code in pressio.
The TPLs will be needed at linking stage when you build *your* code that *uses* pressio.

4. When building your application to use pressio, you just have to point to
the installation directory ``<where-you-cloned-pressio>/install`` with the installed
pressio headers, and you can access all pressio functionalities via the C++ include ``#include<pressio.hpp>``:

.. code-block:: cpp

    #include "pressio/what_you_need.hpp"
    // ...
    int main(){
     // do something
    }

.. warning::

    The procedure above is advised because it enables ``pressio``
    to properly process the cmake options and turn on/off based
    on certain conditions (as explained above).
    The alternative way to use pressio would be to just clone the repo,
    point to the ``<where-you-cloned-pressio>/include`` subdirectory
    and use cmake variables directly when building your code.
    However, this could have unexpected consequences since
    you would be resposible to set the variables correctly but you would not
    know exactly all the possible constraints.
