.. include:: ../mydefs.rst

.. role:: cpp(code)
   :language: cpp

(Affine) Trial Subspace
=======================

Header: ``<pressio/rom_subspaces.hpp>``

API
---

.. code-block:: cpp

  namespace pressio { namespace rom{

  template<class ReducedStateType, class BasisType, class FullStateType>
  /*impl defined*/ create_trial_subspace(BasisType && basis,
                                         FullStateType && offset,
					 bool isAffine);

  }} // end namespace pressio::rom


Templates and Parameters
~~~~~~~~~~~~~~~~~~~~~~~~

- ``ReducedStateType``: the type to represent a state in the subspace

- ``basis``: the subspace's basis object

- ``offset``: full state instance optionally representing the affine shift

- ``isAffine``: boolean value to indicate whether the space is affine or not:

  - if ``isAffine == true``, then ``offset`` is used as shift vector of the affine subspace

  - if ``isAffine == false``, then ``offset`` is copied-from wherever we need to clone it,
    so in this case it is used purely as a "reference/skeleton" to construct a clone

Constraints
~~~~~~~~~~~

- ``ReducedStateType``: must be either an Eigen vector or a Kokkos rank-1 view

- ``BasisType``: must meet the requirements of |link1|_

- ``FullStateType``: must meet the requirements of |link1|_

.. |link1| replace:: *CopyConstructible*
.. _link1: https://en.cppreference.com/w/cpp/named_req/CopyConstructible

Mandates
~~~~~~~~

- ``BasisType`` and ``FullStateType`` must have the same scalar type, i.e.,
  the following must be true: ``std::is_same< typename pressio::Traits<std::decay_t<BasisType>>::scalar_type,
  typename pressio::Traits<std::decay_t<FullStateType>>::scalar_type >::value``


Preconditions
~~~~~~~~~~~~~

- ``basis`` and ``offset`` must be compatible in shape: for example, if ``basis`` is an
  instance of a matrix class representing a column space and ``offset`` is an instance of a vector class,
  then the extent of ``offset`` must be equal to the number of rows of ``basis``

- ``auto basisClone = pressio::ops::clone(basis)`` must *not* have view semantics, i.e.,
  it must create a new, fully independent instance such that modifying ``basisClone``
  does not affect ``basis``

- ``auto offsetClone = pressio::ops::clone(offset)`` must *not* have view semantics, i.e.,
  it must create a new, fully independent instance such that modifying ``offsetClone``
  does not affect ``offset``

Return value
~~~~~~~~~~~~

An instance of an implementation-defined class that represents an abtraction of a trial subspace.
Use ``auto space = create_trial_subspace(/*args*/)`` to let the compiler deduce the type.

Postconditions and Side Effects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- if an rvalue is passed for ``basis`` or ``offset``, the constructor of the return object
  will try to use move semantics. If move semantics are implemented, the temporary is moved from
  and no new memory allocation occurs

- if an lvalue is passed for ``basis`` or ``offset``, when constructing the return object,
  both ``basis`` and ``offset`` will be cloned via ``pressio::ops::clone()``,
  therefore a new allocation is always performed

- the semantics imply that the returned object always *owns* an instance
  of both ``basis`` and ``offset`` that are either the arguments originally
  passed, or identical deep-copies of them

- if ``isAffine == false``, the returned object is guaranteed to
  meet the ``TrialSubspace`` `concept <rom_concepts/c7.html>`__

- if ``isAffine == true``, the returned object is guaranteed to
  meet the ``AffineTrialSubspace`` `concept <rom_concepts/c8.html>`__.


|

Required kernels
----------------

If you are using one of the data types already supported
internally by pressio then you have nothing to do since the
following kernels are already available and visibile to the compiler.
If you use arbitrary data type, then you have to specialize these kernels
for your types, and ensure you include them so that the compiler can find them.

.. code-block:: cpp

  namespace pressio { namespace ops{

  /*FullStateType*/ clone(const /*FullStateType*/ & fomStateIn)
  {
    // create and return a clone of fomStateIn
  }

  /*BasisType*/ clone(const /*BasisType*/ & objIn)
  {
    // create and return a clone of objIn
  }

  template<class ScalarType>
  void fill(/*FullStateType*/ & fullState, ScalarType value)
  {
    // fill the full state with value
  }

  template<class ScalarType>
  void update(/*FullStateType*/ & y, const ScalarType & alpha,
	      const /*FullStateType*/ & x, const ScalarType & beta)
  {
    // compute:
    //   y = alpha*y + beta*x
  }

  template<class AlphaType, class BetaType>
  void product(::pressio::nontranspose /*tag*/,
	       const AlphaType & alpha,
	       const /*BasisType*/ & basis,
	       const /*ReducedStateType*/ & operand,
	       const BetaType & beta,
		/*FullStateType*/ & fullState)
  {
    // compute:
    //   fullState = beta*fullState + alpha * basis * operand
  }

  }} // end namespace pressio::ops
