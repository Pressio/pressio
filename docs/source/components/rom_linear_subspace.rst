.. include:: ../mydefs.rst

.. role:: cpp(code)
   :language: cpp

``LinearSubspace``
==================

Header: ``<pressio/rom_subspaces.hpp>``

Namespace: ``namespace pressio::rom``

|

.. cpp:class:: template <class ReducedStateType, class BasisType, class FullStateType> LinearSubspace

Class template abstracting a linear (optionally affine) subspace.

:tparam ReducedStateType: reduced state's type.
			  Must be either an Eigen vector or a Kokkos rank-1 view.

			  *Requires:* :cpp:`pressio::is_vector_eigen<ReducedStateType>::value == true`
			  *or* :cpp:`pressio::is_vector_kokkos<ReducedStateType>::value == true`

:tparam BasisType: basis type

		   *Requires:* :cpp:`std::is_copy_constructible<BasisType>::value == true`


:tparam FullStateType: full state type.

		       *Requires:* :cpp:`std::is_copy_constructible<FullStateType>::value == true`


*Mandates*

``std::is_same<
typename pressio::Traits<basis_type>::scalar_type,
typename pressio::Traits<full_state_type>::scalar_type>::value == true``

|

.. rubric:: Nested Type Aliases

.. cpp:type:: reduced_state_type

   Same as ``ReducedStateType``.

.. cpp:type:: basis_type

   The ``BasisType`` stripped of any const/reference qualification,
   i.e. ``basis_type = pressio::mpl::remove_cvref_t<BasisType>``

.. cpp:type:: full_state_type

   The ``FullStateType`` stripped of any const/reference qualification,
   i.e. ``full_state_type = pressio::mpl::remove_cvref_t<FullStateType>``

..
   .. rubric:: Private Data Members

   .. cpp:member:: const basis_type LinearSubspace::basis_

      The basis object.

   .. cpp:member:: const full_state_type LinearSubspace::offset_

      The offset object.

   .. cpp:member:: const bool LinearSubspace::isAffine_

      The flag to signal if the subspace is affine or not.


|

.. rubric:: Constructors

.. cpp:function:: LinearSubspace(const basis_type & basis, \
		  const full_state_type & offset, \
		  bool isAffine)

  Constructor accepting lvalues arguments and a boolean.

  :param basis: basis object
  :param offset: affine offset
  :param isAffine: boolean to signal an affine space

  *Preconditions*:

  - ``auto basisClone = pressio::ops::clone(basis)`` must *not* have view semantics,
    i.e., the return value ``basisClone`` must be a new, independent object such that
    modifying it does not affect on the operand object ``basis``

  - ``auto offsetClone = pressio::ops::clone(offset)`` must *not* have view semantics,
    i.e., the return value ``offsetClone`` must be a new, independent object such that
    modifying it does not affect on the operand object ``offset``

  *Semantics and Side Effects*:

  performs new allocations to hold identical deep-copies of the
  arguments ``basis`` and ``offset`` by calling ``pressio::ops::clone()`` for both.

  *Post-conditions*:

  the constructed object owns *immutable* data members which are clones of the arguments
  passed to the constructor


.. cpp:function:: LinearSubspace(basis_type && basis, \
		  full_state_type && offset, \
		  bool isAffine)

  Constructor accepting rvalues arguments and a boolean.

  :param basis: basis object
  :param offset: affine offset
  :param isAffine: boolean to signal an affine space

  *Semantics and Side Effects*:

  The constructor will try to use move semantics. If ``basis_type`` and ``full_state_type``
  implement move semantics, the temporaries arguments are moved from and no new memory allocation occurs.
  If those types do not implement move semantics, then the move operation falls back to a copy,
  and a new allocation occurs.

  *Post-conditions*:

  the constructed object owns *immutable* data members that are identical to the arguments
  passed to the constructor


.. cpp:function:: LinearSubspace(const basis_type & basis, \
		  full_state_type && offset, \
		  bool isAffine)

  Constructor accepting an lvalue and rvalue arguments for basis and offset, and a boolean.

.. cpp:function:: LinearSubspace(basis_type && basis, \
		  const full_state_type & offset, \
		  bool isAffine)

  Constructor accepting an rvalue and lvalue arguments for basis and offset, and a boolean.

|

.. rubric:: Member functions

.. cpp:function:: reduced_state_type createReducedState() const

   Create an instance of a reduce state.

.. cpp:function:: full_state_type createFullState() const

   Create an instance of a full state.

.. cpp:function:: void mapFromReducedState(const reduced_state_type & from, full_state_type & to) const

   Apply the linear transformation to map ``from`` to the corresponding full state ``to``.

.. cpp:function::  full_state_type createFullStateFromReducedState(const reduced_state_type & from) const

   Create a full state instance and map ``from`` to it.
   This is equivalent to doing:

   .. code-block:: cpp

      auto fomStateTo = createFullState();
      mapFromReducedState(from, fomStateTo);


.. cpp:function:: const basis_type & viewBasis() const

   Vew the basis.
