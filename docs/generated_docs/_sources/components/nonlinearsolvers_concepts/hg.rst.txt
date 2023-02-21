.. role:: raw-html-m2r(raw)
   :format: html

.. include:: ../../mydefs.rst

``SystemWithHessianAndGradient``
================================

.. literalinclude:: ../../../../include/pressio/solvers_nonlinear/concepts/solvers_system_hessian_gradient.hpp
   :language: cpp
   :lines: 52-101

Semantic requirements
---------------------

:red:`finish`


Syntax-only snippet
-------------------

.. code-block:: cpp

   struct Something
   {
     using state_type    = /* your type */;
     using hessian_type  = /* your type */;
     using gradient_type = /* your type */;
     using residual_norm_type = /* your type */;

     state_type    createState() const;
     hessian_type  createHessian() const;
     gradient_type createGradient() const;

     void residualNorm(const state_type & state,
                       pressio::Norm normKind,
                       residual_norm_type & resNorm) const;

     void gradient(const state_type &,
                   gradient_type &,
                   pressio::Norm normKind,
                   residual_norm_type & normResidual,
                   bool recomputeJacobian) const;

     void hessian(const state_type &, hessian_type &) const;
   };


 ..
    //
    // execution/language axioms
    //
    axiom NonAliasingOperators(const T & A){
     auto s1 = A.createState();
     auto s2 = A.createState();
     std::addressof(s1) != std::addressof(s2);

     auto g1 = A.createGradient();
     auto g2 = A.createGradient();
     std::addressof(g1) != std::addressof(g2);

     auto H1 = A.createHessian();
     auto H2 = A.createHessian();
     std::addressof(H1) != std::addressof(H2);
    } &&
    //
    axiom DeterministicResidualNorm(){
      // residual norm is always deterministic
    } &&
    axiom DeterministicHessian(){
      // hessian computation is always deterministic
    } &&
    axiom DeterministicGradientOnRequest(){
      // the C++ standard actually names Deterministic as EqualityPreserving
      // gradient is deterministic (or equalitypreserving) iff recomputeJacobian == true
    } &&
    axiom BlockingOperations(){
      // all methods are blocking (allocations/operations complete before returning)
    } &&
    axiom ConstCorrectness(){
      // const qualification is preserved, methods do NOT modify const arguments
    } &&

    //
    // mathematical axioms
    //
    requires(){
      requires RealVectorSpaceElement<typename T::state_type>;
      requires RealVectorSpaceElement<typename T::gradient_type>;
      requires RealVectorSpaceElement<typename T::hessian_type>;
    } &&
    axiom Overdetermined(){
      // dimension of the state vector space > dimension of the residual vector space
      // i.e. # of equations > # of unknowns
    } &&
    axiom SumOfSquaredResidualsCostFunction(){
      // we are minimizing ||F||^2
    } &&
    axiom FullRankJacobian(){
      // # of rows > # cols, so full rank if the matrix is full column rank
    } &&
    axiom SecondOrderHessianTruncation(){
      // hessian  = J^T J
      // so neglecting sum_j r_j grad^2 r_j
    } &&
    axiom DifferentiableResidual(){
      // residual function is supposed to be differentiable
    };
