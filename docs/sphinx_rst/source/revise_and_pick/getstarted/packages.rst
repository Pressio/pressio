.. role:: raw-html-m2r(raw)
   :format: html


Components
==========

The pressio C++ library is divided into several components:

.. list-table::
   :header-rows: 1

   * - Name
     - Description
     - Link
     - Reference header(s)
   * - @m_span{m-text m-success}mpl@m_endspan
     - metaprogramming functionalities
     - `Source <https://github.com/Pressio/pressio/tree/main/include/mpl>`_
     - ``#include<pressio_mpl.hpp>``
   * - @m_span{m-text m-success}utils@m_endspan
     - common functionalities\ :raw-html-m2r:`<br/>`\ e.g., I/O helpers, static constants, etc
     - `Source <https://github.com/Pressio/pressio/tree/main/include/utils>`_
     - ``#include<pressio_utils.hpp>``
   * - @m_span{m-text m-success}containers@m_endspan
     - wrappers for vectors, matrices and multi-vectors, and expressions (span, diagonal and subspan)
     - `Source <https://github.com/Pressio/pressio/tree/main/include/containers>`_
     - ``#include<pressio_containers.hpp>``
   * - @m_span{m-text m-success}ops@m_endspan
     - shared-memory and distributed linear algebra
     - `Source <https://github.com/Pressio/pressio/tree/main/include/ops>`_
     - ``#include<pressio_ops.hpp>``
   * - @m_span{m-text m-success}apps@m_endspan
     - suites of mini-apps used for basic testing
     - `Source <https://github.com/Pressio/pressio/tree/main/include/apps>`_
     - ``#include<pressio_apps.hpp>``
   * - @m_span{m-text m-success}qr@m_endspan
     - QR factorization functionalities
     - `Source <https://github.com/Pressio/pressio/tree/main/include/qr>`_
     - ``#include<pressio_qr.hpp>``
   * - @m_span{m-text m-success}solvers@m_endspan
     - linear and non-linear solvers :raw-html-m2r:`<br>` (e.g., Newton-Raphson, Gauss-Newton, Levenberg-Marquardt)
     - `Source <https://github.com/Pressio/pressio/tree/main/include/solvers>`_
     - ``#include<pressio_solvers.hpp>``
   * - @m_span{m-text m-success}ode@m_endspan
     - explicit only methods :raw-html-m2r:`<br/>`\ implict only methods :raw-html-m2r:`<br/>` all
     - :raw-html-m2r:`<br/>`\ :raw-html-m2r:`<br/>`\ `Source <https://github.com/Pressio/pressio/tree/main/include/ode>`_
     - ``#include<pressio_ode_explicit.hpp>``\ :raw-html-m2r:`<br/>` ``#include<pressio_ode_implicit.hpp>`` :raw-html-m2r:`<br/>` ``#include<pressio_ode.hpp>``
   * - @m_span{m-text m-success}rom@m_endspan
     - Galerkin ROMs :raw-html-m2r:`<br/>` LSPG ROMs :raw-html-m2r:`<br/>` WLS ROMs :raw-html-m2r:`<br/>` all
     - :raw-html-m2r:`<br/>`\ :raw-html-m2r:`<br/>`\ :raw-html-m2r:`<br/>`\ `Source <https://github.com/Pressio/pressio/tree/main/include/rom>`_
     - ``#include<pressio_rom_galerkin.hpp>`` :raw-html-m2r:`<br/>` ``#include<pressio_rom_lspg.hpp>`` :raw-html-m2r:`<br/>` ``#include<pressio_rom_wls.hpp>`` :raw-html-m2r:`<br/>` ``#include<pressio_rom.hpp>``


The top-down order used above is informative of the dependency structure.
For example, every package depends on ``mpl``. The ``ops`` package depends only on ``mpl``\ , ``utils``\ , ``containers``.
At the bottom of the stack we have the ``rom`` package which requires all the others.

This structure of the framework has several benefits.


* 
  Maintability: ``pressio`` can be more easily developed and maintained since its components depend on one another through well-defined public interfaces,
  and appropriate namespaces are used to organize classes.

* 
  Selective usability: this modular framework allows users, if needed, to leverage invidual components.
  For instance, if a user needs/wants just the QR methods, they can simply use that package,
  and all the dependencies on the others are enabled automatically.

@m_class{m-block m-warning}

@par
When you use functionalities from a specific package, you should just include
the corresponding header and the dependencies (based on the explanation above) are included automatically.
For example, if you want to do Galerkin with explicit time integration,
you just do ``#include <pressio_rom_galerkin.hpp>`` because all the needed
packages are automatically included. There is not need to manually include all of them yourself.
In the future, we might refine further the granularity of the headers to allow a finer control.

@m_class{m-block m-warning}

@par One header to include them all
If you want to access *all* functionalities, you can use:

.. code-block:: cpp

   #include "pressio.hpp"
