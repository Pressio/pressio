.. include:: ../../mydefs.rst

``ComposableIntoDefaultProblem``
================================

Header: ``<pressio/rom_concepts.hpp>``

.. important::

   To avoid name conflicts, this belongs to the namespace: ``pressio::rom::galerkin::steady``

.. literalinclude:: ../../../../include/pressio/rom/concepts/galerkin_steady_default.hpp
   :language: cpp
   :lines: 54-98

The ``DefaultProblem`` concept specifies that an object
of type ``TrialSubspaceType`` and an object of type ``FomSystemType`` can
be composed to define a steady default Galerkin problem.


Semantic requirements
---------------------

:red:`finish`

..
   The concept is modeled only if it is satisfied,
   all subsumed concepts are modeled and letting
   ``S`` be an instance of ``SubspaceType``,
   then all of the following hold:

   - a reduced residual instance ``rGal`` has extent equal to the dimensionality of the subspace, i.e.
     ``pressio::ops::extent(rGal, 0) == S.dimension()``

   - a reduced dense jacobian instance ``JGal`` is such that
     ``pressio::ops::extent(JGal, 0) == pressio::ops::extent(JGal, 1) == S.dimension()``
     must be true, i.e. ``JGal`` must be a square matrix of size
     equal to the dimensionality of the subspace

   - *blocking operations*: all methods are blocking

   - *equality preserving*: equal inputs yield equal outputs

   - *const correctness*: methods may modify only the non-constant operands.
     Operands that are constant must not be modified.
