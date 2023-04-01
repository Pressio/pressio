.. role:: raw-html-m2r(raw)
   :format: html

Gauss-Newton (via QR)
=====================

Header: ``<pressio/solvers_nonlinear_gaussnewton.hpp>``

API
---

.. literalinclude:: ../../../include/pressio/solvers_nonlinear/solvers_create_gauss_newton.hpp
   :language: cpp
   :lines: 54, 232-233, 245-266, 268, 313,316


Parameters
~~~~~~~~~~

.. list-table::
   :widths: 18 82
   :header-rows: 1
   :align: left

   * -
     -

   * - ``system``
     - your problem instance

   * - ``qrSolver``
     - qr solver to use within each nonlinear iteration

Constraints
~~~~~~~~~~~

With C++20, the constraints would be enforced via concepts using
the *requires-clause* shown in the API synopsis above.
Since we cannot yet officially upgrade to C++20, the constraints
are currently enforced via static asserts (to provide a decent error message) and/or SFINAE.

The concepts are documented `here <nonlinsolvers_concepts.html>`__.

Examples
--------

TBD
