.. include:: ../mydefs.rst

.. role:: cpp(code)
   :language: cpp

``TrialColumnSubspace``
=======================

Header: ``<pressio/rom_subspaces.hpp>``

Namespace: ``namespace pressio::rom``

API
---

.. code-block:: cpp

   namespace pressio{ namespace rom{

   template <class ReducedStateType, class BasisType, class FullStateType>
   class TrialColumnSubspace
   {
     using reduced_state_type = ReducedStateType;
     using basis_type         = std::remove_cv_t<BasisType>;
     using full_state_type    = std::remove_cv_t<FullStateType>;
     using offset_type        = full_state_type;

     TrialColumnSubspace(const basis_type & basis,                                     (1)
			 const full_state_type & offset,
			 bool isAffine);

     TrialColumnSubspace(basis_type && basis,                                          (2)
			 full_state_type && offset,
			 bool isAffine);

     TrialColumnSubspace(const basis_type & basis,                                     (3)
			 full_state_type && offset,
			 bool isAffine);

     TrialColumnSubspace(basis_type && basis,                                          (4)
			 const full_state_type & offset,
			 bool isAffine);

     reduced_state_type createReducedState() const;                                    (5)

     full_state_type createFullState() const;                                          (6)

     void mapFromReducedState(const ReducedStateType &,                                (7)
		              full_state_type &) const;

     full_state_type createFullStateFromReducedState(const ReducedStateType &) const;  (8)

     const basis_type & viewBasis() const;                                             (9)

     const offset_type & viewAffineOffset() const;                                    (10)

     bool isAffine() const;                                                           (11)
   };

   }} //end namespace

Constraints
-----------

- :cpp:`pressio::is_vector_eigen<ReducedStateType>::value == true ||
  pressio::is_vector_kokkos<ReducedStateType>::value == true`

- :cpp:`std::is_copy_constructible<basis_type>::value == true &&
  std::is_pointer<basis_type>::value == false &&
  pressio::mpl::is_std_shared_ptr<basis_type>::value == false`

- :cpp:`std::is_copy_constructible<full_state_type>::value == true &&
  std::is_pointer<full_state_type>::value == false &&
  pressio::mpl::is_std_shared_ptr<full_state_type>::value == false`

Mandates
--------

- :cpp:`std::is_same<typename pressio::Traits<basis_type>::scalar_type,
  typename pressio::Traits<full_state_type>::scalar_type>::value == true`


Preconditions
-------------

- ``basis`` *represents a column space*. Note that here we could have used
  a more specific name like ``basisMatrix`` but we intentionally kept the
  a more generic name becuase we want to allow users to fully exploit the power of templates here.
  Doing so allows us to support various scenarios: the most basic one is when you
  provide ``basis`` to be an instance of a ``basis_type`` class representing
  a "concrete" matrix (e.g. Eigen::Matrix, Tpetra::MultiVector, etc),
  with its columns containing the basis vectors. A more flexible scenario is one where
  ``basis`` is an instance of a class that is NOT a "concrete matrix" but only an
  abstraction. This is useful for example if you did not have a way to store the matrix,
  but you can only implement the action of the matrix.

- for 1,3, the following expression ``auto basisClone = pressio::ops::clone(basis)``
  must *not* have view semantics, i.e., ``basisClone`` must be an object independent
  of ``basis`` such that any modifications applied to ``basisClone`` do not affect
  the original ``basis`` object

- for 1,4, the following expression ``auto offsetClone = pressio::ops::clone(offset)``
  must *not* have view semantics, i.e., ``offsetClone`` must be an object independent
  of ``offset`` such that any modifications applied to ``offsetClone`` do not affect
  the original ``offset`` object

- if you pass ``true`` to the ``isAffine`` constructor parameter, then ``offset``
  is expected to be non trivial offset. If you pass true and ``offset`` is trivial,
  i.e. represents a "zero" offset, then it is as if you wanted a non-affine subspace
  but are using an affine implementation

Semantics and Side Effects
---------------------------

- constructor 1: makes new allocations by calling ``pressio::ops::clone(...)`` on
  both ``basis`` and ``offset``, and stores the returned objects as data members

- constructor 2: this will try to use move semantics. If ``basis_type`` and ``full_state_type``
  implement move semantics, the temporary arguments are moved from and
  no new memory allocation occurs. If those types do not implement move semantics,
  then the move operation falls back to a copy, and a new allocation occurs.
  The returned objects will own these new objects as members.

- constructors 3,4: mixed behavior for lvalues and rvalues arguments, see 1,2.


Postconditions
--------------

- the constructed object owns *immutable* data members which are clones (according to ``ops::clone()``)
  of the arguments passed to the constructor


Non-member functions
--------------------

The following free function is provided as an alternative
to directly instantiating the class:

.. code-block:: cpp

   namespace pressio{ namespace rom{

   template<class ReducedStateType, class BasisType, class FullStateType>
   auto create_trial_column_subspace(BasisType && basis,
                                     FullStateType && offset,
				     bool isAffine)

   }} //end namespace
