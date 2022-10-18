.. role:: raw-html-m2r(raw)
   :format: html

Gauss-Newton via normal-equations
=================================

Defined in header ``<pressio/solvers_nonlinear.hpp>``

API
---

.. code-block:: cpp

   namespace pressio{ namespace nonlinearsolvers{

   template<class SystemType, class LinearSolverType>
   #ifdef PRESSIO_ENABLE_CXX20
     requires
	  (OverdeterminedSystemWithResidualAndJacobian<SystemType>
	|| OverdeterminedSystemWithFusedResidualAndJacobian<SystemType>
	|| SystemWithHessianAndGradient<SystemType>
	|| SystemWithFusedHessianAndGradient<SystemType>)
       && LinearSolverForNonlinearLeastSquares<
	    mpl::remove_cvref_t<LinearSolverType>,
	    typename SystemType::state_type>
   #endif
   auto create_gauss_newton(const SystemType & system,          (1)
                            LinearSolverType && lsolver);

   template<
     class SystemType,
     class LinearSolverType,
     class WeightingOpType>
   #ifdef PRESSIO_ENABLE_CXX20
     requires
	  (OverdeterminedSystemWithResidualAndJacobian<SystemType>
	|| OverdeterminedSystemWithFusedResidualAndJacobian<SystemType>)
       && LinearSolverForNonlinearLeastSquares<
	    mpl::remove_cvref_t<LinearSolverType>,
	    typename SystemType::state_type>
       && LeastSquaresWeightingOperator<
	    mpl::remove_cvref_t<WeightingOpType>,
	    typename SystemType::residual_type,
	    typename SystemType::jacobian_type>
   #endif
   auto create_gauss_newton(const SystemType & system,          (2)
                            LinearSolverType && lsolver,
                            WeightingOpType && weightOperator);

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

   * - ``weightOperator``
     - the weighting operator for doing weighted least-squares.

Constraints
~~~~~~~~~~~

Each overload is associated with a set of constraints.
With C++20, these would be enforced via concepts using
the *requires-clause* shown in the API synopsis above.
Since we cannot yet officially upgrade to C++20, the constraints
are currently enforced via static asserts (to provide a decent error message)
and/or SFINAE. The concepts used are:

- `nonlinearsolvers::OverDeterminedSystemWithResidualAndJacobian <nonlinearsolvers_concepts/rj_ovdet.html>`__

- `nonlinearsolvers::OverDeterminedSystemWithFusedResidualAndJacobian <nonlinearsolvers_concepts/rj_fused_ovdet.html>`__

- `nonlinearsolvers::SystemWithHessianAndGradient <nonlinearsolvers_concepts/hg.html>`__

- `nonlinearsolvers::SystemWithFusedHessianAndGradient <nonlinearsolvers_concepts/hg_fused.html>`__

- `nonlinearsolvers::LinearSolverForNonlinearLeastSquares <nonlinearsolvers_concepts/c4.html>`__


Example usage
^^^^^^^^^^^^^

bla blas
