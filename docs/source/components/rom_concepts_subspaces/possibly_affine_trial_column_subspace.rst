
.. role:: raw-html(raw)
    :format: html

``PossiblyAffineTrialColumnSubspace``
=====================================

Header: ``<pressio/rom_concepts.hpp>``

Namespace: ``pressio::rom``

.. literalinclude:: ../../../../include/pressio/rom/concepts/possibly_affine_trial_column_subspace.hpp
   :language: cpp
   :lines: 56-95


Semantic requirements
---------------------

The concept is modeled only if it is satisfied,
all subsumed concepts are modeled and given
an instance ``s``, of type ``T``, all of the following hold:

- *non aliasing instantiation*: given the following:

  .. code-block:: cpp

     auto r1 = s.createReducedState();
     auto r2 = s.createReducedState();

  ``r1`` and ``r2`` must be distinct objects, ``std::addressof(r1) != std::addressof(r2)``,
  and such that any modification to ``r1`` does not affect ``r2``, and the same applies
  to ``createFullState``.

- *blocking operations*: all methods are blocking, meaning that all temporary
  allocations and operations are complete before the methods return
  and not outstanding work remains

- *const correctness*: methods may modify only the non-constant operands.
  Operands that are constant must not be modified.

- ``s.isColumnSpace() == true`` and ``s.isRowSpace() == false`` always

- ``auto dim = s.dimension()`` represents the true dimensionality of the subspace, :raw-html:`<br/>`
  i.e. ``dim == pressio::ops::extent(s.basisOfTranslatedSpace(), 1)`` is true

- if the subspace *is affine*, then:

  - ``auto & shift = s.translationVector()`` should return a non-zero
    translation identifying the non trivial origin of the affine subspace

  - ``auto & basis = s.basisOfTranslatedSpace()`` returns the basis of the translated space

  - calling ``s.basis()`` is undefined behavior

- if the subspace *is NOT affine*, then:

  - ``auto & shift = s.translationVector()`` returns a trivial translation,
    namely a translation with all "zeros" representing a zero offset

  - ``auto & basis = s.basis()`` and ``auto & basis2 = s.basisOfTranslatedSpace()``
    are equivalent


..
   Syntax only
   -----------

   .. literalinclude:: ./syntax_only_subspaces_concepts.cc
      :language: cpp
      :lines: 18-35



..
      template <class T, class Enable = void>
      struct CanOverwriteFullStateArgument : std::true_type{};

      template <class T>
      struct CanOverwriteFullStateArgument<
	T,
	std::enable_if_t<
	  std::is_void<
	    decltype(
	     std::declval<T const>().mapFromReducedState(int{}, std::declval< int&& >())
	     )
	    >::value
	  >
	> : std::false_type{};

      } end namespace impl
