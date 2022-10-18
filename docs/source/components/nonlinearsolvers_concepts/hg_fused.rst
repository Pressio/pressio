.. role:: raw-html-m2r(raw)
   :format: html

.. include:: ../../mydefs.rst

``SystemWithFusedHessianAndGradient``
=====================================

.. literalinclude:: ../../../../include/pressio/solvers_nonlinear/concepts/solvers_system_fused_hessian_gradient.hpp
   :language: cpp
   :lines: 52-99

Semantic requirements
---------------------

:red:`finish`


Syntax-only snippet
-------------------

.. code-block:: cpp

   struct ProblemHessGradFused
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

     void hessianAndGradient(const state_type &,
			     hessian_type &,
			     gradient_type &,
			     pressio::Norm normKind,
			     residual_norm_type & normResidual,
			     bool recomputeJacobian) const;
   };
