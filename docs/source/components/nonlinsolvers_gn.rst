.. role:: raw-html-m2r(raw)
   :format: html

Gauss-Newton
============

.. note::

    Defined in header ``<pressio/solvers_nonlinear.hpp>``

    Public namespace: ``pressio::nonlinearsolvers``

Gauss-Newton via Normal-Equations
---------------------------------

API, Parameters and Requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   template<
     class ProblemClassType,
     class StateType,
     class LinearSolverType
     >                                                             (1)
   auto create_gauss_newton(const ProblemClassType & system,
                            const StateType & state,
                            LinearSolverType && lsolver);

   template<
     class ProblemClassType,
     class StateType,
     class LinearSolverType,
     class WeightingOpType
     >                                                             (2)
   auto create_gauss_newton(const ProblemClassType & system,
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

----

Gauss-Newton via QR factorization
---------------------------------

API, Parameters and Requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   template<
     class ProblemClassType,
     class StateType,
     class QRSolverType
     >
   auto create_gauss_newton(const ProblemClassType & system,
                            const StateType & state,
                            QRSolverType && qrsolver);

* 
  ``system``\ :

  * 
    instance of the problem you want to solve

    .. warning::

        * *only* accepts the `residual-jacobian API or its fused version <nonlinsolvers_system_api.html>`_

* 
  ``state``\ :

  * your state data
  * Requirements: must be an Eigen or Kokkos vector: \todo explain why

* 
  ``qrsolver``\ :

  * solver needed to solve the QR-based formulation of the least-squares problem `see this <https://en.wikipedia.org/wiki/QR_decomposition>`_
  * we suggest to use the `pressio QR package <qr.html>`_
  * if you want to implement your own, then it has to conform to the `this API <qr.html>`_

Ops
^^^

Example usage
^^^^^^^^^^^^^
