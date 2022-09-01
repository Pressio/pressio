.. role:: raw-html-m2r(raw)
   :format: html

.. include:: mydefs.rst

pressio C++ library
===================

.. admonition:: :medium:`Advancing reduced order models (ROMs) for dynamical systems in science and engineering`
    :class: note

    This is the documentation of the `C++ library <https://github.com/Pressio/pressio>`__, one element of the `Pressio ecosystem <https://pressio.github.io/>`_.

.. :summarylineindexpage:`Advancing reduced order models (ROMs) for dynamical systems in science and engineering.`

.. This is the documentation of the `C++ library <https://github.com/Pressio/pressio>`__, one component of the `Pressio ecosystem <https://pressio.github.io/>`_.


The following table represents the *software stack* of the functionalities in pressio:
each component (level) in the stack depends on all the ones below it.
This structure is mirrored in the API documentation (see the left sidebar).

.. list-table::
   :widths: 10 48 42
   :header-rows: 1
   :align: left

   * -
     - Description
     - Header(s)

   * - ``rom``
     - linear subspaces :raw-html-m2r:`<br/>` Galerkin: steady :raw-html-m2r:`<br/>` Galerkin: unsteady :raw-html-m2r:`<br/>` LSPG: steady :raw-html-m2r:`<br/>` LSPG: unsteady :raw-html-m2r:`<br/>` :raw-html-m2r:`<br/>`
     - ``<pressio/rom_subspaces.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/rom_galerkin_steady.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/rom_galerkin_unsteady.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/rom_lspg_steady.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/rom_lspg_unsteady.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/rom.hpp>`` :small:`to include all`

   * - ``ode``
     - explicit steppers :raw-html-m2r:`<br/>` implicit steppers :raw-html-m2r:`<br/>` ``advance_<keywords>`` functions :raw-html-m2r:`<br/>` :raw-html-m2r:`<br/>`
     - ``<pressio/ode_steppers_explicit.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/ode_steppers_implicit.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/ode_advancers.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/ode.hpp>`` :small:`to include all`

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

   * - ``concepts``
     - common concepts
     - ``<pressio/concepts.hpp>``

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

* `overview <overview.html>`_: motivations, objectives and how have built it from scratch

* `how to install <installation.html>`_: it is a header-only library, should be trivial

* `explore the tutorials <https://pressio.github.io/pressio-tutorials>`_


Generic programming and concepts
--------------------------------

Arguably the main foundation of pressio is the use
of generic programming--*or, more humbly, we can
at least say that it is what we strive for*.
Since the early development stages, we have tried to rely on
concept-driven design. Here, the term concept does not necessarily
refer to the C++ concepts feature introduced in C++20.
Rather, you can think of it more loosely as "what properties/syntax a type meets,
what you can do with it and, also, what a type should definitely satisfy".
So, at least for now, think of the term concept a bit loosely.
In pressio, in fact, we cannot yet upgrade to C++20, and thus
we cannot use and *stricto sensu* enforce concepts.
For now, concepts are present but enforced via, e.g., SFINAE,
or in some cases (abusing the syntax) with static asserts.
Later on, when feasible, we will enforce them as the standard dictates.

The message we want to convey is that *"concepts" are a fundamental
design part of pressio and are present everywhere*---note, that, if you have
used C++ templates, you *have* used concepts knowingly or not.
In our documentation, we make the effort to highlight the use of concepts
by dedicating to each component of the library a full section
(see the sidebar menu on the left)
to discuss and formalize how concepts are used in that component.




License and Citation
--------------------

The full license (BSD-3) is available `here <https://pressio.github.io/various/license/>`_.

We are working on publishing this: you can find our arXiv preprint at: https://arxiv.org/abs/2003.07798

Questions?
----------

Find us on Slack: https://pressioteam.slack.com or
open an issue on `github <https://github.com/Pressio/pressio>`_.


.. toctree::
   :maxdepth: 1
   :hidden:

   overview
   installation

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
   ./components/concepts
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
