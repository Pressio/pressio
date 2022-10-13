.. role:: raw-html-m2r(raw)
   :format: html

Gauss-Newton (via normal-equations)
===================================

Defined in header ``<pressio/solvers_nonlinear.hpp>``


API, Parameters and Requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   namespace pressio{ namespace nonlinearsolvers{

   template<class ProblemClassType, class LinearSolverType>
   auto create_gauss_newton(const ProblemClassType & system,    (1)
                            LinearSolverType && lsolver);

   template<class ProblemClassType, class LinearSolverType, class WeightingOpType>
   auto create_gauss_newton(const ProblemClassType & system,    (2)
                            LinearSolverType && lsolver,
                            WeightingOpType && weightOperator);

   }}

* ``system``

  - instance of your problem class defining the problem

    .. warning::

        * overload 1 accepts a system satisfying *any of* the following concepts:
	  `SystemWithResidualAndJacobian <nonlinearsolvers_concepts/c1.html>`__,
	  `SystemWithFusedResidualAndJacobian <nonlinearsolvers_concepts/c2.html>`__,
	  `SystemWithHessianAndGradient <nonlinearsolvers_concepts/c3.html>`__
	  or `SystemWithFusedHessianAndGradient <nonlinearsolvers_concepts/c3.html>`__

        * overload 2 *only* accepts a system satisfying either the
	  `SystemWithResidualAndJacobian <nonlinearsolvers_concepts/c1.html>`__ or
	  `SystemWithFusedResidualAndJacobian <nonlinearsolvers_concepts/c2.html>`__

* ``lsolver``

  * linear solver for solving the normal equations, choose one from `linear solvers <linsolvers.html>`_
  * if you want to implement your own, then it has to conform to the `linear solver API <linsolvers.html>`_

* ``weightOperator``:

  * weighting operator for doing weighted least-squares.
    Must conform to:

    .. code-block:: cpp

       class WeightingOperator
       {
         public:
         void operator()(const residual_type & operand, residual_type & result);
         void operator()(const jacobian_type & operand, jacobian_type & result);
       };


Example usage
^^^^^^^^^^^^^

bla blas
