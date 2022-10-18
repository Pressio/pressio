Newton-Raphson
==============

Defined in header ``<pressio/solvers_nonlinear.hpp>``

API
---

.. code-block:: cpp

   namespace pressio{ namespace nonlinearsolvers{

   template<class SystemType, class LinearSolverType>
   #ifdef PRESSIO_ENABLE_CXX20
     requires
	(DeterminedSystemWithResidualAndJacobian<SystemType>
      || DeterminedSystemWithFusedResidualAndJacobian<SystemType>)
       && LinearSolverForNewtonRaphson<
	    mpl::remove_cvref_t<LinearSolverType>,
	    typename SystemType::jacobian_type,
	    typename SystemType::residual_type,
	    typename SystemType::state_type>
   #endif
   auto create_newton_raphson(const SystemType & system,
                              LinearSolverType && lsolver);

   }}

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

   * - ``lsolver``
     - linear solver to use within each nonlinear iteration


Constraints
~~~~~~~~~~~

Each overload is associated with a set of constraints.
With C++20, these would be enforced via concepts using
the *requires-clause* shown in the API synopsis above.
Since we cannot yet officially upgrade to C++20, the constraints
are currently enforced via static asserts (to provide a decent error message)
and/or SFINAE. The concepts used are:

- `nonlinearsolvers::DeterminedSystemWithResidualAndJacobian <nonlinearsolvers_concepts/rj.html>`__

- `nonlinearsolvers::DeterminedSystemWithFusedResidualAndJacobian <nonlinearsolvers_concepts/rj_fused.html>`__

- `nonlinearsolvers::LinearSolverForNewtonRaphson <nonlinearsolvers_concepts/c4.html>`__


Example usage
-------------

.. code-block:: cpp

   int main()
   {
     // assuming that:
     // problem_t is a problem class that meets API
     // state_t is defined too

     problem_t myProblem;

     // create linear system
     using lin_solver_t = /* something that meets API needed */;
     lin_solver_t linearSolverObj;

     namespace pnls = pressio::nonlinearsolvers;
     auto NonLinSolver = pnls::create_newton_raphson(myProblem, linearSolverObj);

     state_t y(10);
     // set initial state
     NonLinSolver.solve(myProblem, y);
   }
