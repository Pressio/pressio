
Levenberg-Marquardt
===================

Header: ``<pressio/solvers_nonlinear_levmarq.hpp>``

API
---

.. literalinclude:: ../../../include/pressio/solvers_nonlinear/solvers_create_levenberg_marquardt.hpp
   :language: cpp
   :lines: 54-55, 69-89, 143-144

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

   * - ``linSolver``
     - linear solver to use within each nonlinear iteration

   * - ``weigher``
     - the weighting operator for doing weighted least-squares

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
