.. role:: raw-html-m2r(raw)
   :format: html

.. include:: ../../mydefs.rst

``DeterminedSystemWithResidualAndJacobian``
===========================================

.. literalinclude:: ../../../../include/pressio/solvers_nonlinear/concepts/solvers_system_residual_jacobian.hpp
   :language: cpp
   :lines: 93-100

Semantic requirements
---------------------

:red:`finish`

- dimension of the residual vector space == dimension of the state vector space
  i.e. # of equations is same as # of unknowns
