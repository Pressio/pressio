.. include:: ../../mydefs.rst

``SemiDiscreteFom``
===================

Header: ``<pressio/rom_concepts.hpp>``

.. literalinclude:: ../../../../include/pressio/rom/concepts/fom_semi_discrete.hpp
   :language: cpp
   :lines: 58-100

Semantic requirements
---------------------

Given an instance ``A`` of type ``T``, ``SemiDiscreteFom<T>``
is modeled if it is satisfied, all subsumed concepts are modeled and:

- all methods are blocking, meaning that all temporary
  allocations and operations needed to execute those methods
  are completed and no outstanding work remains upon return

- methods may modify only the non-constant operands.
  Operands that are constant must not be modified.

- ``auto rhs = A.createRightHandSide()``

  - returns an object with all its "elements" zero initialized

- doing:

  .. code-block:: cpp

     auto rhs1 = A.createRightHandSide();
     auto rhs2 = A.createRightHandSide();
     //...
     auto rhsN = A.createRightHandSide();

  implies that ``rhs1, rhs2, ..., rhsN`` must be distinct objects,
  and such that any modification to ``rhs1`` does not affect any of the others
  and viceversa.
  In other words, calling ``A.createRightHandSide()`` yields independent instances.

- ``A.rightHandSide(state, evalTime, rhs)``

  - overwrites ``rhs`` with the result

  - is equality preserving, i.e. given equal
    inputs ``state, evalTime``, the result written to ``rhs`` remains the same
