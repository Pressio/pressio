.. role:: raw-html-m2r(raw)
   :format: html

.. include:: ../../mydefs.rst

``OverdeterminedSystemWithFusedResidualAndJacobian``
====================================================

.. literalinclude:: ../../../../include/pressio/solvers_nonlinear/concepts/solvers_system_fused_residual_jacobian.hpp
   :language: cpp
   :lines: 105-112

Semantic requirements
---------------------

:red:`finish`

- dimension of the residual vector space > dimension of the state vector space
  i.e. # of equations > # of unknowns
