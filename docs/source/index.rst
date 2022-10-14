.. role:: raw-html-m2r(raw)
   :format: html

.. include:: mydefs.rst

pressio C++ library
===================

.. admonition:: :medium:`Advancing reduced order models for dynamical systems in science and engineering`
    :class: note

    This is the documentation of the `C++ library <https://github.com/Pressio/pressio>`__, one element of the `Pressio ecosystem <https://pressio.github.io/>`_.


This work was started with a focus on projection-based reduced-order models (ROMs),
which is a strongly multidisciplinary topic.
Working towards a production-level ROM capability inevitably means touching
multiple fields, ranging from, e.g., linear algebra, nonlinear solvers
and optimization, to time integration, distributed computing and HPC.
The complexity level is further amplified if the goal is to develop a *generic* library.
**Modularity, abstractions and well-defined APIs** thus become fundamental design principles
needed to properly build such a project from the ground up.

This has been part of our effort from the beginning, and has lead to the following
"stacked" design of pressio: each component (level) handles a specific capability and depends,
via well-defined public APIs, on the ones below it. This approach has several benefits, but
the main one is that each component becomes usable on its own, and, as a whole,
the stack constitutes the foundation of the top-level ``pressio/rom`` component.
The API documentation mirros this structure.

..
  The following table represents the *software stack* of the functionalities in pressio:
  each component (level) in the stack depends on all the ones below it.


.. list-table::
   :widths: 10 48 42
   :header-rows: 1
   :align: left

   * -
     - Description
     - Header(s)

   * - ``rom``
     - (linear) subspaces :raw-html-m2r:`<br/>` Galerkin: steady :raw-html-m2r:`<br/>` Galerkin: unsteady :raw-html-m2r:`<br/>` LSPG: steady :raw-html-m2r:`<br/>` LSPG: unsteady :raw-html-m2r:`<br/>` :raw-html-m2r:`<br/>`
     - ``<pressio/rom_subspaces.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/rom_galerkin_steady.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/rom_galerkin_unsteady.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/rom_lspg_steady.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/rom_lspg_unsteady.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/rom.hpp>`` :small:`includes all above`

   * - ``ode``
     - explicit steppers :raw-html-m2r:`<br/>` implicit steppers :raw-html-m2r:`<br/>` ``advance_<keywords>`` functions :raw-html-m2r:`<br/>` :raw-html-m2r:`<br/>`
     - ``<pressio/ode_steppers_explicit.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/ode_steppers_implicit.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/ode_advancers.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/ode.hpp>`` :small:`includes all above`

   * - ``solvers_nonlinear``
     - e.g., Newton-Raphson, Gauss-Newton, Levenberg-Marquardt.
     - ``<pressio/solvers_nonlinear.hpp>``

   * - ``solvers_linear``
     - linear dense (on-node) solvers
     - ``<pressio/solvers_linear.hpp>``

   * - ``qr``
     - QR factorization functionalities
     - ``<pressio/qr.hpp>``

   * - ``ops``
     - shared-memory/distributed linear algebra kernels specializations
     - ``<pressio/ops.hpp>``

   * - ``expressions``
     - expressions templates, e.g.: span, diagonal, subspan
     - ``<pressio/expressions.hpp>``

   * - ``type_traits``
     - type traits and detection
     - ``<pressio/type_traits.hpp>``

   * - ``utils``
     - logger, constants, etc
     - ``<pressio/utils.hpp>``

   * - ``mpl``
     - metaprogramming functionalities
     - ``<pressio/mpl.hpp>``


Get Started
-----------

* `how to install <installation.html>`_: it is a header-only library, should be trivial

* `explore the tutorials <https://pressio.github.io/pressio-tutorials>`_


Generic programming and concepts
--------------------------------

Arguably the main foundation of pressio is the use of
generic programming--*or, more humbly, we can at least say that it is what we strive for*.
Since the early development stages, we have relied on concept-driven design.
Here, the term concept does not necessarily
refer to the C++ concepts feature introduced in C++20.
You can think of it more broadly as "what properties/syntax a type meets,
what you can do with it and, also, what a type should definitely satisfy".

Note, that, if you have used or use C++ templates, you *have* used
concepts knowingly or not. This is becuase when you write a function or class
template, you have some expectations of what a template is going to do.
Concepts are basically a way to *explicitly* formalize those expectations.

..
   The message we want to convey is that *"concepts" are a fundamental
   design part of pressio*. In our documentation, we make the effort to
   highlight the use of concepts
   by dedicating to each component of the library a full section
   to discuss and formalize how concepts are used in that component.

Until we can stably upgrade to C++20, we cannot properly use C++20 concepts.
The concepts in pressio are guarded inside a preprocessor directive.
To enable them, it suffices to use a C++20 compiler and set ``-DCMAKE_CXX_STANDARS=20``.
When disabled, the pressio code by default enforces them via, e.g., SFINAE,
or in some cases (abusing the syntax) with static asserts.



License and Citation
--------------------

The full license (BSD-3) is available `here <https://github.com/Pressio/pressio/blob/main/LICENSE>`_.

We are working on publishing this: you can find our arXiv preprint at: https://arxiv.org/abs/2003.07798

Questions?
----------

Find us on Slack: https://pressioteam.slack.com or
open an issue on `github <https://github.com/Pressio/pressio>`_.


.. toctree::
   :maxdepth: 1
   :hidden:

   installation
   keywords

.. toctree::
   :caption: API
   :maxdepth: 1
   :hidden:

   ./components/rom
   ./components/ode
   ./components/nonlinsolvers
   ./components/linsolvers
   ./components/qr
   ./components/ops
   ./components/expressions
   ./components/type_traits
   ./components/utils
   ./components/mpl

.. toctree::
   :caption: Miscellanea
   :maxdepth: 1
   :hidden:

   Tutorials <https://pressio.github.io/pressio-tutorials>
   GitHub Repo <https://github.com/Pressio/pressio>
   Open an issue/feature req. <https://github.com/Pressio/pressio/issues>
   license
