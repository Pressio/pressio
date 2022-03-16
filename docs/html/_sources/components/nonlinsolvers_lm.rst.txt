Levenberg-Marquardt
===================

.. note::

    Defined in header ``<pressio/solvers_nonlinear.hpp>``

    Public namespace: ``pressio::nonlinearsolvers``

Levenberg–Marquardt
-------------------

API, Parameters and Requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   template<
     class ProblemClassType,
     class StateType,
     class LinearSolverType
     >                                                             (1)
   auto create_levenberg_marquardt(const ProblemClassType & system,
                                   const StateType & state,
                                   LinearSolverType && lsolver);

   template<
     class ProblemClassType,
     class StateType,
     class LinearSolverType,
     class WeightingOpType
     >                                                             (2)
   auto create_levenberg_marquardt(const ProblemClassType & system,
                                   const StateType & state,
                                   LinearSolverType && lsolver,
                                   WeightingOpType && weightOperator);

* 
  ``system``\ :

  * 
    instance of the problem you want to solve

    .. warning::

        * overload 1: accepts the `residual-jacobian, hessian-gradient API, or their fused versions <nonlinsolvers_system_api.html>`_
        * overload 2: *only* accepts the `residual-jacobian API or its fused version <nonlinsolvers_system_api.html>`_

* 
  ``state``\ :

  * your state data
  * Requirements: must be an Eigen or Kokkos vector: \todo explain why

* 
  ``lsolver``\ :

  * linear solver for solving the normal equations, choose one from `linear solver API <linsolvers.html>`_
  * if you want to implement your own, then the linear solver class still has to conform to the `linear solver API <linsolvers.html>`_

* 
  ``weightOperator``\ :

  * weighting operator for doing weighted least-squares.
    Must conform to:

    .. code-block:: cpp

       class WeightingOperator
       {
       public:
       void operator()(const residual_type & operand, residual_type & result);
       void operator()(const jacobian_type & operand, jacobian_type & result);
       };

Ops
^^^

Example usage
^^^^^^^^^^^^^
