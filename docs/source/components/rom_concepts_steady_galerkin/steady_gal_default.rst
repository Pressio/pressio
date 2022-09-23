
``ProjectableOnPossiblyAffineSubspace``
=======================================

Header: ``<pressio/rom_concepts.hpp>``

.. important::

   This concept belongs to the namespace: ``pressio::rom::galerkin::steady``

.. literalinclude:: ../../../../include/pressio/rom/concepts/steady_galerkin_projectable_on_affine_subspace.hpp
   :language: cpp
   :lines: 59-106

Semantic requirements
---------------------

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
