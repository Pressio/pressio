.. role:: raw-html-m2r(raw)
   :format: html

Iteratively reweighted least squares
====================================

Defined in header ``<pressio/solvers_nonlinear.hpp>``


API
---

.. code-block:: cpp

   namespace pressio{ namespace nonlinearsolvers{ namespace experimental{

   template<class ProblemClassType, class LinearSolverType>
   #ifdef PRESSIO_ENABLE_CXX20
     requires
	  (OverdeterminedSystemWithResidualAndJacobian<SystemType>
	|| OverdeterminedSystemWithFusedResidualAndJacobian<SystemType>)
       && LinearSolverForNonlinearLeastSquares<
	    mpl::remove_cvref_t<LinearSolverType>,
	    typename SystemType::state_type>
   #endif
   auto create_irls_gauss_newton(const ProblemClassType & system,
		                 LinearSolverType && lsolver);

   }}}


.. warning::

   Note that this is still inside the experimental namespace.


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

- `nonlinearsolvers::OverDeterminedSystemWithResidualAndJacobian <nonlinearsolvers_concepts/rj_ovdet.html>`__

- `nonlinearsolvers::OverDeterminedSystemWithFusedResidualAndJacobian <nonlinearsolvers_concepts/rj_fused_ovdet.html>`__

- `nonlinearsolvers::LinearSolverForNonlinearLeastSquares <nonlinearsolvers_concepts/c4.html>`__


Example usage
^^^^^^^^^^^^^

something
