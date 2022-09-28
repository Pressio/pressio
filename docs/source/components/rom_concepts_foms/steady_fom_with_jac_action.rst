.. include:: ../../mydefs.rst

``SteadyFomWithJacobianAction``
===============================

Header: ``<pressio/rom_concepts.hpp>``

.. literalinclude:: ../../../../include/pressio/rom/concepts/fom_steady_with_jac_action.hpp
   :language: cpp
   :lines: 58-86

Semantic requirements
---------------------

Given an instance ``A`` of type ``T`` and an object ``operand``
of type ``JacobianActionOperandType``,
``SteadyFomWithJacobianAction<T, JacobianActionOperandType>``
is modeled if it is satisfied, all subsumed concepts are modeled and:

- all methods are blocking, meaning that all temporary
  allocations and operations needed to execute those methods
  are completed and no outstanding work remains upon return

- methods may modify only the non-constant operands.
  Operands that are constant must not be modified.

- ``auto r = A.createResidual()`` and
  ``auto result = A.createApplyJacobianResult(operand)`` return objects
  with all "elements" zero initialized

- doing:

  .. code-block:: cpp

     auto r1 = A.createResidual();
     auto r2 = A.createResidual();
     //...
     auto rN = A.createResidual();

  implies that ``r1, r2, ..., rN`` must be distinct objects,
  and such that any modification to ``r1`` does not affect any of the others
  and viceversa.
  In other words, calling ``A.createResidual()`` yields independent instances.
  And similarly applies to ``A.createApplyJacobianResult()``.

- ``A.residual(state, r)`` and ``A.applyJacobian(state, operand, ja)``

  - overwrite ``r`` and ``ja`` with their respective results

  - both are equality preserving, i.e. equal inputs imply equal outputs

  - let ``J`` represent the Jacobian which we compute the action of,
    then ``J`` must be the jacobian of the residual evaluated for the
    same ``state``. In other words, the Jacobian used for
    computing its action must be mathematically "consistent" with the residual.
