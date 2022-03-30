.. role:: raw-html-m2r(raw)
   :format: html

.. include:: mydefs.rst

pressio C++ library
===================

.. admonition:: Advancing reduced order models (ROMs) for dynamical systems in science and engineering
    :class: note

    This is the documentation of the `C++ library <https://github.com/Pressio/pressio>`__, one component of the `Pressio ecosystem <https://pressio.github.io/>`_.

.. :summarylineindexpage:`Advancing reduced order models (ROMs) for dynamical systems in science and engineering.`

.. This is the documentation of the `C++ library <https://github.com/Pressio/pressio>`__, one component of the `Pressio ecosystem <https://pressio.github.io/>`_.


.. list-table::
   :widths: 10 50 40
   :header-rows: 1
   :align: left

   * -
     - Description
     - Header(s)

   * - :packnameindexpage:`mpl`
     - metaprogramming functionalities
     - ``<pressio/mpl.hpp>``

   * - :packnameindexpage:`utils`
     - logger, constants, etc
     - ``<pressio/utils.hpp>``

   * - :packnameindexpage:`type_traits`
     - type traits and detection
     - ``<pressio/type_traits.hpp>``

   * - :packnameindexpage:`expressions`
     - Expr templates for span, diagonal, subspan, etc.
     - ``<pressio/expressions.hpp>``

   * - :packnameindexpage:`ops`
     - specializations of shared-memory and distributed linear algebra kernels
     - ``<pressio/ops.hpp>``

   * - :packnameindexpage:`qr`
     - QR factorization functionalities
     - ``<pressio/qr.hpp>``

   * - :packnameindexpage:`solvers_linear`
     - linear *dense* solvers
     - ``<pressio/solvers_linear.hpp>``

   * - :packnameindexpage:`solvers_nonlinear`
     - Including Newton-Raphson, Gauss-Newton, Levenberg-Marquardt
     - ``<pressio/solvers_nonlinear.hpp>``

   * - :packnameindexpage:`ode`
     - explicit time steppers :raw-html-m2r:`<br/>` implicit time steppers :raw-html-m2r:`<br/>` advancers :raw-html-m2r:`<br/>` :raw-html-m2r:`<br/>`
     - ``<pressio/ode_steppers_explicit.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/ode_steppers_implicit.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/ode_advancers.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/ode.hpp>`` :small:`to include all`

   * - :packnameindexpage:`rom`
     - decoder :raw-html-m2r:`<br/>` Galerkin :raw-html-m2r:`<br/>` LSPG: steady/unsteady :raw-html-m2r:`<br/>` WLS :raw-html-m2r:`<br/>` :raw-html-m2r:`<br/>`
     - ``<pressio/rom_decoder.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/rom_galerkin.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/rom_lspg.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/rom_wls.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/rom.hpp>`` :small:`to include all`


.. raw:: html

   <!-- This structure benefits: -->
   <!-- * Maintability: the components depend on one another through well-defined public interfaces, -->
   <!-- and appropriate namespaces are used to properly scope the code. -->

.. raw:: html

   <!-- * Selective usability: users can leverage invidual functionalities. -->
   <!-- This allows finer-grained control on what you include and use. -->

Get Started
-----------

*
  `introduction <introduction.html>`_: describes objectives and design ideas

*
  `how to install <installation.html>`_: it is a header-only library, should be trivial

*
  `explore the tutorials <https://pressio.github.io/pressio-tutorials>`_

.. raw:: html

   <!-- ## What if your types are not natively supported in pressio? -->

.. raw:: html

   <!-- Check if your types are supported by lookig at the -->
   <!-- [dependencies](md_pages_getstarted_build_and_install.html): if they are -->
   <!-- listed there, most likely you are good to go, and you don't need to provide extra information to pressio. -->

.. raw:: html

   <!-- Not supported? You can file an [issue](https://github.com/Pressio/pressio/issues) -->
   <!-- to request it and wait on it, or can proceed -->
   <!-- as in [tutorialsB](./md_pages_tutorials_tutorial1udops.html). Or do both! -->

License and Citation
--------------------

The full license (BSD-3) is available `here <https://pressio.github.io/various/license/>`_.

We are working on publishing this: you can find our arXiv preprint at: https://arxiv.org/abs/2003.07798

Questions?
----------

Find us on Slack: https://pressioteam.slack.com or
open an issue on `github <https://github.com/Pressio/pressio>`_.


|


.. toctree::
   :maxdepth: 1

   introduction
   installation
   GitHub Repo <https://github.com/Pressio/pressio>
   Open an issue/feature req. <https://github.com/Pressio/pressio/issues>
   license
   ./components/mpl
   ./components/utils
   ./components/type_traits
   ./components/expressions
   ./components/ops
   ./components/qr
   ./components/linsolvers
   ./components/nonlinsolvers
   ./components/ode
   ./components/rom
