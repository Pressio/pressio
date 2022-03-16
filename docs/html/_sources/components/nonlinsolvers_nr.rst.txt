Newton-Raphson
==============

.. note::

    Defined in header ``<pressio/solvers_nonlinear.hpp>``

    Public namespace: ``pressio::nonlinearsolvers``

API, Parameters and Requirements
--------------------------------

.. code-block:: cpp

   template<class ProblemClassType, class StateType, class LinearSolverType>
   auto create_newton_raphson(const ProblemClassType & system,
                              const StateType & state,
                              LinearSolverType && lsolver);

* 
  ``system``\ :

  * instance of the problem you want to solve
  * Requirements: must conform to either the residual-jacobian API or the fused residual-jacobian API

* 
  ``state``\ :

  * your state
  * Requirements: must be copy constructible

* 
  ``lsolver``\ :

  * linear solver called at each inner iteration
  * Requirements: must conform to this `API <linsolvers.html>`_

Ops
---

When using custom data types not supported in `pressio ops <ops.html>`_\ , you need to specialize a trait class and some operations
such that pressio can operate on your data. For the sake of explanation, suppose that you use:

.. code-block:: cpp

   using scalar_type   = double;
   using state_type    = ACustomStateType;    //this can be any type
   using jacobian_type = ACustomJacobianType; //this can be any type

Then you need to provide the following specializations:

.. code-block:: cpp

   namespace pressio{

   template<> struct Traits<ACustomStateType>{
     using scalar_type = double;
   };

   template<> struct Traits<ACustomJacobianType>{
     using scalar_type = double;
   };

   namespace ops{

   std::size_t extent(ACustomStateType &, int i){
     /* return extent along i-th dimension */
   }

   std::size_t extent(ACustomJacobianType &, int i){
     /* return extent along i-th dimension */
   }

   ACustomStateType    clone(const ACustomStateType & src){
     /* return a deep copy of src */
   }

   ACustomJacobianType clone(const ACustomJacobianType & src){
     /* return a deep copy of src */
   }

   void set_zero(ACustomStateType & object){    /* set elements zero */ }
   void set_zero(ACustomJacobianType & object){ /* set elements zero */ }

   scalar_type norm2(const ACustomStateType & object){
     /* return l2-norm of object */
   }

   void update(ACustomStateType & v,        scalar_type a,
               const ACustomStateType & v1, scalar_type b)
   {
     // compute v = a*v + v1*b;
   }

   void scale(ACustomStateType & v, scalar_type factor){
     /* scale v elementwise by factor
   }

   }}//end namepsace pressio::ops

Example usage
-------------

.. code-block:: cpp

   int main()
   {
     // assuming that:
     // problem_t is a problem class that meets API
     // state_t is defined too

     problem_t myProblem;
     state_t y(10);
     // set initial state

     // create linear system
     using lin_solver_t = /* something that meets API needed */;
     lin_solver_t linearSolverObj;

     namespace pnls = pressio::nonlinearsolvers;
     auto NonLinSolver = pnls::create_newton_raphson(myProblem, y, linearSolverObj);
     NonLinSolver.solve(myProblem, y);
   }
