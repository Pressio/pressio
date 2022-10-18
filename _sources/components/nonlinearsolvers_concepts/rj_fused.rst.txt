.. role:: raw-html-m2r(raw)
   :format: html

.. include:: ../../mydefs.rst

``SystemWithFusedResidualAndJacobian``
======================================

.. literalinclude:: ../../../../include/pressio/solvers_nonlinear/concepts/solvers_system_fused_residual_jacobian.hpp
   :language: cpp
   :lines: 52-93

Semantic requirements
---------------------

:red:`finish`

Syntax-only snippet
-------------------

.. code-block:: cpp

   struct ProblemFusedResJac
   {
     using state_type     = /* your type */;
     using residual_type  = /* your type */;
     using jacobian_type  = /* your type */;

     state_type    createState() const;
     residual_type createResidual() const;
     jacobian_type createJacobian() const;
     void residualAndJacobian(const state_type& x,
                              residual_type & res,
                              jacobian_type & jac,
                              bool recomputeJacobian) const;
   };
