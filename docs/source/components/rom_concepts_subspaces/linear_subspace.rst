
``VectorSubspace``
==================

Header: ``<pressio/rom_concepts.hpp>``

Namespace: ``pressio::rom``

.. literalinclude:: ../../../../include/pressio/rom/concepts/linear_subspace.hpp
   :language: cpp
   :lines: 56-67

Semantic requirements
---------------------

The concept is modeled only if it is satisfied,
all subsumed concepts are modeled and given
an instance ``s``, of type ``T``, all of the following hold:

- ``s`` is immutable and, consequently, so are the underlying basis
  and subspace it represents.

- let ``auto & basis = s.basis()``, then:

  - if ``s.isColumnSpace() == true``, then ``basis`` is full *column* rank

  - if ``s.isRowSpace() == true``, then ``basis`` is full *row* rank

- if ``s.isColumnSpace() == true``, then it ``s.isRowSpace() == false`` and vice versa

- ``auto dim = s.dimension()`` represents the true dimensionality of the subspace

- must be closed under addition and scalar multiplication

..
   Syntax only
   -----------

   .. literalinclude:: ./syntax_only_subspaces_concepts.cc
      :language: cpp
      :lines: 8-16
