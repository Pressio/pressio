.. include:: ../../mydefs.rst

``SemiDiscreteFom``
===================

Header: ``<pressio/rom_concepts.hpp>``

Namespace: ``pressio::rom``

.. literalinclude:: ../../../../include/pressio/rom/concepts/semi_discrete_fom.hpp
   :language: cpp
   :lines: 56-75

Semantic requirements
---------------------

:red:`finish`

..
   - *non aliasing instantiation*: given the following:

     .. code-block:: cpp

	auto r1 = A.createRightHandSide();
	auto r2 = A.createRightHandSide();

     ``r1`` and ``r2`` must be distinct objects, ``std::addressof(r1) != std::addressof(r2)``,
     and such that any modification to ``r1`` does not affect ``r2``

   - *blocking operations*: all methods are blocking, meaning that all temporary
     allocations and operations are complete before the methods return and not outstanding work remains

   - *equality preserving*: given ``A`` an object of type `T`, calling ``A.rightHandSide(...)``
     with equal inputs yields equal outputs.

   - *const correctness*: methods may modify only the non-constant operands.
     Operands that are constant must not be modified.


.. Syntax only
.. -----------

.. .. literalinclude:: ./syntax_only_fom_system_concepts.cc
..    :language: cpp
..    :lines: 6-18
