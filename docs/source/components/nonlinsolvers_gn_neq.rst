.. role:: raw-html-m2r(raw)
   :format: html

Gauss-Newton (via normal eqs)
=============================

Defined in ``<pressio/solvers_nonlinear_gaussnewton.hpp>``

API
---

.. literalinclude:: ../../../include/pressio/solvers_nonlinear/solvers_create_gauss_newton.hpp
   :language: cpp
   :lines: 54-55, 65-84, 53,53,53, 140-167, 303-304

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
are currently enforced via static asserts (to provide a decent error message)
and/or SFINAE. The concepts are documented `here <nonlinearsolvers_concepts>`__.

Return Value
~~~~~~~~~~~~

Returns a solver object exposing the following public API:

.. code-block:: cpp

   // This is not the actual class, it just describes the API
   class Solver
   {

   public:
     // query/set update criterion
     Update currentUpdateCriterion() const;
     void setUpdateCriterion(Update value);

     // query/set stop criterion, tolerance
     Stop currentStopCriterion() const;
     void setStopCriterion(Stop value);
     void setStopTolerance(ScalarType value);
     void setMaxIterations(int newMax);

     void solve(StateType & solutionInOut);
   };

Examples
^^^^^^^^

.. admonition:: Demos
   :class: tip

   1. `full demo for the Rosenbrock function <https://pressio.github.io/pressio-tutorials/using_eigen/nonlinsolvers1.html>`__
