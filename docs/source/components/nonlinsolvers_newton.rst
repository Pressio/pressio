Newton-Raphson
==============

Header: ``<pressio/solvers_nonlinear_newton.hpp>``

API
---

.. literalinclude:: ../../../include/pressio/solvers_nonlinear/solvers_create_newton.hpp
   :language: cpp
   :lines: 59-60, 78-108, 152-153

Parameters
~~~~~~~~~~

.. list-table::
   :widths: 18 82
   :header-rows: 1
   :align: left

   * -
     -

   * - ``system``
     - instance of your problem

   * - ``linSolver``
     - linear solver to use within each nonlinear iteration


Constraints
~~~~~~~~~~~

With C++20, the constraints would be enforced via concepts using
the *requires-clause* shown in the API synopsis above.
Since we cannot yet officially upgrade to C++20, the constraints
are currently enforced via static asserts (to provide a decent error message) and/or SFINAE.

The concepts are documented `here <nonlinsolvers_concepts.html>`__.

..
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
--------

TBD: link to tutorial example


..
   - `nonlinearsolvers::DeterminedSystemWithResidualAndJacobian <nonlinearsolvers_concepts/rj.html>`__

   - `nonlinearsolvers::DeterminedSystemWithFusedResidualAndJacobian <nonlinearsolvers_concepts/rj_fused.html>`__

   - `nonlinearsolvers::LinearSolverForNewtonRaphson <nonlinearsolvers_concepts/c4.html>`__
