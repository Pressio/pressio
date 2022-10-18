.. role:: raw-html-m2r(raw)
   :format: html

Gauss-Newton via QR factorization
===================================

Defined in header ``<pressio/solvers_nonlinear.hpp>``


API, Parameters and Requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   namespace pressio{ namespace nonlinearsolvers{

   template<class SystemType, class QRSolverType>
   #ifdef PRESSIO_ENABLE_CXX20
     requires
	  (OverdeterminedSystemWithResidualAndJacobian<SystemType>
	|| OverdeterminedSystemWithFusedResidualAndJacobian<SystemType>)
       && QRSolverForGnQr<
	    mpl::remove_cvref_t<LinearSolverType>,
	    typename SystemType::state_type,
	    typename SystemType::jacobian_type,
	    typename SystemType::residual_type>
   #endif
   auto create_gauss_newtonQR(const SystemType & system,
                              QRSolverType && qrsolver);

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

   * - ``qrsolver``
     - solver to use within each nonlinear iteration

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

- `nonlinearsolvers::QRSolverForGnQr <nonlinearsolvers_concepts/c4.html>`__

Example usage
^^^^^^^^^^^^^

bla
