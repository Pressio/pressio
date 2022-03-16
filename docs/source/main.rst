.. role:: raw-html-m2r(raw)
   :format: html

Pressio C++ Library
===================

.. note::
    *Advancing reduced order models (ROMs) for dynamical systems in science and engineering.*

    This is the documentation of the `C++ library <https://github.com/Pressio/pressio>`_\ , one component of the `Pressio ecosystem <https://pressio.github.io/>`_.

.. list-table::
   :header-rows: 1

   * - Name
     - Description/Content
     - Links
     - Corresponding header(s)
   * - mpl
     - :raw-html-m2r:`<br/>` metaprogramming functionalities
     - `Code <https://github.com/Pressio/pressio/tree/develop/include/pressio/mpl>`_ :raw-html-m2r:`<br/>` `Documentation <components/mpl.html>`_
     - ``<pressio/mpl.hpp>``
   * - utils
     - :raw-html-m2r:`<br/>` logging, constants, helpers, etc
     - `Code <https://github.com/Pressio/pressio/tree/develop/include/pressio/utils>`_\ :raw-html-m2r:`<br/>`\ `Documentation <components/utils.html>`_
     - ``<pressio/utils.hpp>``
   * - type_traits
     - :raw-html-m2r:`<br/>` traits/detection classes
     - `Code <https://github.com/Pressio/pressio/tree/develop/include/pressio/type_traits>`_\ :raw-html-m2r:`<br/>`\ `Documentation <components/type_traits.html>`_
     - ``<pressio/type_traits.hpp>``
   * - expressions
     - :raw-html-m2r:`<br/>` classes for various abstractions (span, diagonal, subspan, etc.)
     - `Code <https://github.com/Pressio/pressio/tree/develop/include/pressio/expressions>`_\ :raw-html-m2r:`<br/>`\ `Documentation <components/expressions.html>`_
     - ``<pressio/expressions.hpp>``
   * - ops
     - :raw-html-m2r:`<br/>` specializations of shared-memory and distributed linear algebra kernels
     - `Code <https://github.com/Pressio/pressio/tree/develop/include/pressio/ops>`_\ :raw-html-m2r:`<br/>`\ `Documentation <components/ops.html>`_
     - ``<pressio/ops.hpp>``
   * - qr
     - :raw-html-m2r:`<br/>` QR factorization functionalities
     - `Code <https://github.com/Pressio/pressio/tree/develop/include/pressio/qr>`_\ :raw-html-m2r:`<br/>`\ `Documentation <components/qr.html>`_
     - ``<pressio/qr.hpp>``
   * - solvers_linear
     - :raw-html-m2r:`<br/>` linear solvers (wrappers around existing TPLs)
     - `Code <https://github.com/Pressio/pressio/tree/develop/include/pressio/solvers_linear>`_\ :raw-html-m2r:`<br/>`\ `Documentation <components/linsolvers.html>`_
     - ``<pressio/solvers_linear.hpp>``
   * - solvers_nonlinear
     - :raw-html-m2r:`<br/>` general info :raw-html-m2r:`<br/>` Newton-Raphson :raw-html-m2r:`<br/>` Gauss-Newton :raw-html-m2r:`<br/>` Levenberg-Marquardt :raw-html-m2r:`<br/>`
     - `Code <https://github.com/Pressio/pressio/tree/develop/include/pressio/solvers_nonlinear>`_ :raw-html-m2r:`<br/>` `Documentation <components/nonlinsolvers_general.html>`_ :raw-html-m2r:`<br/>` `Documentation <components/nonlinsolvers_nr.html>`_ :raw-html-m2r:`<br/>` `Documentation <components/nonlinsolvers_gn.html>`_ :raw-html-m2r:`<br/>` `Documentation <components/nonlinsolvers_lm.html>`_
     - ``<pressio/solvers_nonlinear.hpp>``
   * - ode
     - :raw-html-m2r:`<br/>` explicit steppers :raw-html-m2r:`<br/>`\ implicit steppers :raw-html-m2r:`<br/>` advancers :raw-html-m2r:`<br/>`
     - `Code <https://github.com/Pressio/pressio/tree/develop/include/pressio>`_ :raw-html-m2r:`<br/>` `Documentation <components/ode_steppers_explicit.html>`_\ :raw-html-m2r:`<br/>` `Documentation <components/ode_steppers_implicit.html>`_ :raw-html-m2r:`<br/>`\ `Documentation <components/ode_advance.html>`_
     - :raw-html-m2r:`<br/>` ``<pressio/ode_steppers_explicit.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/ode_steppers_implicit.hpp>``\ :raw-html-m2r:`<br/>` ``<pressio/ode_advancers.hpp>`` :raw-html-m2r:`<br/>` (for everything: ``<pressio/ode.hpp>``\ )
   * - rom
     - :raw-html-m2r:`<br/>`\ general info :raw-html-m2r:`<br/>` decoder :raw-html-m2r:`<br/>` Galerkin\ :raw-html-m2r:`<br/>` LSPG: steady\ :raw-html-m2r:`<br/>` LSPG: unsteady\ :raw-html-m2r:`<br/>` WLS\ :raw-html-m2r:`<br/>`
     - `Code <https://github.com/Pressio/pressio/tree/develop/include/pressio/rom>`_ :raw-html-m2r:`<br/>`\ `Documentation <components/rom_general.html>`_ :raw-html-m2r:`<br/>`\ `Documentation <components/rom_decoder.html>`_ :raw-html-m2r:`<br/>` `Documentation <components/rom_galerkin.html>`_ :raw-html-m2r:`<br/>` `Documentation <components/rom_lspg_steady.html>`_ :raw-html-m2r:`<br/>` `Documentation <components/rom_lspg_unsteady.html>`_ :raw-html-m2r:`<br/>`  `Documentation <components/rom_wls.html>`_ :raw-html-m2r:`<br/>`
     - :raw-html-m2r:`<br/>` :raw-html-m2r:`<br/>` ``<pressio/rom_decoder.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/rom_galerkin.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/rom_lspg.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/rom_lspg.hpp>``\ :raw-html-m2r:`<br/>` ``<pressio/rom_wls.hpp>`` :raw-html-m2r:`<br/>` (for everything: ``<pressio/rom.hpp>``\ )

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
  `read the introduction <introduction.html>`_ providing an overview, objectives and design ideas

* 
  `how to install <installation.html>`_\ : it is a header-only library, should be trivial

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

.. raw:: html

   <!--
   @m_class{m-note m-success}

   Pressio is an open-source project aimed at enabling leading-edge projection-based
   reduced order models (\proms) for dynamical systems in science and engineering.

   ## Motivation
   Projection-based model reduction refers to a class of surrogate models
   that reduce the number of degrees
   of freedom in the full-order model (FOM) through a projection process.
   This projection step applied to the governing equations often enables one
   to make stronger performance guarantees
   (e.g., of structure preservation, of accuracy via adaptivity) than other
   surrogates like data-fits and perform more accurate *a posteriori*
   error analysis (e.g., via *a posteriori* error bounds or error models).

   Despite these benefits, the practical challenges of
   implementing model-reduction techniques in large-scale codes often
   precludes their adoption in practice; this occurs because standard implementations
   require modifying low-level operations and solvers for each simulation code of interest.
   This implementation strategy is not practical or sustainable
   in many modern settings, because industrial simulation codes often evolve rapidly,
   institutions may employ dozens of simulation codes for different analyses,
   and commercial codes typically do not expose the required low-level
   operators and solvers.


   @m_class{m-note m-success}

   Pressio aims to mitigate the implementation burden of projection-based model
   reduction in large-scale applications without compromising performance.


   ## Main steps of pROMs
   Projection-based model reduction can be broken into three main steps,
   namely data collection, basis creation, and ROM deployment.

   - data collection: \todo (all)

   - compute basis: \todo (all)

   - create/run the ROM: \todo (all)


   @m_class{m-block m-warning}

   @par
   pressioproj currently contains capabilities to perform the last step.
   \todo Say that we have plans for the other steps too.
   Maybe at some point we will provide tools to run the samples,
   but for now that is not a huge priority. we can develop something
   later on to aid this step. For example interfacing with efficient
   POD libraries, providing tools for specific mesh formats (exodus).
    -->

.. raw:: html

   <!--
   ## The Pressio framework
   \pressioproj is a computational *framework*, comprising a (growing) collection of repositories :

   * [pressio](https://github.com/Pressio/pressio): &emsp;&ensp;&emsp;&emsp;&ensp;core C++ library based on generic programming;

   to support applications with arbitrary data types; -->
   <!-- [pressio4py](https://github.com/Pressio/pressio4py): &emsp;&emsp;&nbsp;&nbsp;Python bindings for the core Pressio C++ functionalities; -->
   <!-- [pressio-builder](https://github.com/Pressio/pressio-builder): &nbsp;&nbsp;&nbsp;auxiliary bash scripts for building/testing; -->
   <!-- [pressio-tutorials](https://github.com/Pressio/pressio-tutorials): &nbsp;tutorials explaining how to use `pressio` and its functionalities.

   ## Where to go from here
   If you are new and want to learn more, start from the [userguide](./md_pages_get_started.html)
   and see how to install and use pressio, or you can jump directly
   to the [tutorials](./md_pages_tutorials.html)
   and/or [examples](md_pages_examples.html) -->
