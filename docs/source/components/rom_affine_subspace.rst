.. include:: ../mydefs.rst

.. role:: cpp(code)
   :language: cpp

``AffineSubspace``
========================

Header: ``<pressio/rom_subspaces.hpp>``

API
---

.. code-block:: cpp

   namespace pressio{ namespace rom{

   template <class BasisType, class OffsetType>
   class AffineSubspace
   {
     public:
       // type aliases
       using basis_type  = std::remove_cv_t<BasisType>;
       using offset_type = std::remove_cv_t<OffsetType>;

       // constructors
       AffineSubspace(const basis_type & basis,                             (1)
		      const offset_type & offset,
		      bool isAffine);

       AffineSubspace(basis_type && basis,                                  (2)
		      offset_type && offset,
		      bool isAffine);

       AffineSubspace(const basis_type & basis,                             (3)
		      offset_type && offset,
		      bool isAffine);

       AffineSubspace(basis_type && basis,                                  (4)
		      const offset_type & offset,
		      bool isAffine);

       // copy constructor and copy assignment
       AffineSubspace(const AffineSubspace &);                        (5)
       AffineSubspace & operator=(const AffineSubspace&) = delete;    (6)

       // accessors
       const basis_type & viewBasis() const;                                      (7)

       const offset_type & viewAffineOffset() const;                              (8)

       bool isAffine() const;                                                     (9)
   };

   }} //end namespace

Constraints
-----------

- :cpp:`std::is_copy_constructible<basis_type>::value == true`

- :cpp:`std::is_copy_constructible<offset_type>::value == true`

- :cpp:`std::is_pointer<basis_type>::value == false`

- :cpp:`std::is_pointer<offset_type>::value == false`

- :cpp:`pressio::mpl::is_std_shared_ptr<basis_type>::value == false`

- :cpp:`pressio::mpl::is_std_shared_ptr<offset_type>::value == false`

Mandates
--------

- :cpp:`std::is_same<typename pressio::Traits<basis_type>::scalar_type,
  typename pressio::Traits<offset_type>::scalar_type>::value == true`


Member functions
----------------

Constructors
^^^^^^^^^^^^

:memberfunction:`(1)`: constructs an instance of this class from lvalue arguments

  *Preconditions*:

  - ``basis`` must be *an object of a class representing a linear basis*.
    Note that here we could have limited the meaning,
    e.g. narrowing it to ``basisMatrix``, but we intentionally keep the scope more generic
    becuase it allows users to fully exploit the power of templates.
    Doing so, in fact, allows us to support various scenarios: the most basic one
    is when ``basis`` is an instance of a ``basis_type`` class representing
    a "concrete" matrix (e.g. Eigen::Matrix, Tpetra::MultiVector, etc),
    with its columns or rows containing the basis vectors spanning the space.
    An alternative more flexible (but also most delicate) scenario is one where
    ``basis`` is an instance of a class that is NOT a "concrete matrix" but only
    an abstraction/wrapper of one.

  - the following expression ``auto basisClone = pressio::ops::clone(basis)`` is
    valid with requirements/semantics as `explained here <ops/clone.html>`__, i.e.
    ``basisClone`` is an object independent of ``basis``
    and any operation/modification applied to ``basisClone`` (or any data it references/uses)
    do not affect the original ``basis`` object (or any data it references/uses)

  - the following expression ``auto offsetClone = pressio::ops::clone(offset)`` is
    valid with requirements/semantics as `explained here <ops/clone.html>`__, i.e.
    ``offsetClone`` is an object independent of ``offset``
    and any operation/modification applied to ``offsetClone`` (or any data it references/uses)
    do not affect the original ``offset`` object (or any data it references/uses)

  *Effects*:

  - the constructor performs the following operations:
    it calls ``auto ret = pressio::ops::clone(...)`` on both ``basis``
    and ``offset``, and uses these returned object to copy-construct
    const-qualified data members, therefore new allocations
    are made to store the data members

  *Postconditions*:

  - the constructed object is immutable and owns, until its destruction,
    const-qualified data members that are clones (according to ``ops::clone()``)
    of the arguments originally passed to the constructor. Therefore, this class
    has a clear invariant that is established upon construction and cannot be
    changed until its destruction.

  - the immutability and invariant of an instance of this class relies
    on the assumption that users do not violate const correctness or
    abuse of const_casting.


:memberfunction:`(2)`: constructs an instance of this class from rvalues arguments

  *Preconditions*:

  - ``basis`` has the same preconditions as in (1).

  *Effects*:

  - this constructor performs the following: it tries to use move semantics
    by move-constructing ``basis`` and ``offset``. If ``basis_type`` and/or
    ``offset_type`` implement move semantics, the temporary arguments are
    moved-constructed and no new memory allocation occurs to allocate the data members.
    If move semantics are not implemented, then the move operation
    falls back to a copy, and a new allocation (potentially) occurs.

  *Postconditions*:

  - the constructed object is immutable and owns, until its destruction,
    const-qualified data members that were either move constructed from the
    original arguments, or are just copies. As in (1), this class
    has a clear invariant that is established upon construction and cannot be
    changed until its destruction.


:memberfunction:`(3)`: constructs an instance of this class from an lvalue and rvalue

:memberfunction:`(4)`: constructs an instance of this class from an rvalue and lvalue


:memberfunction:`(5)`: copy constructor

  - The semantics of the copy constructor are the same as in constructor (1):
    regardless of what ``basis_type`` and ``offset_type`` are, invoking the copy
    constructor *always* creates a new instance whose members are "cloned" from ``other``.

:memberfunction:`(6)`: copy assignment operator is deleted

   - this is explicitly marked as deleted because it would violate
     the class immutability


Accessors
^^^^^^^^^

:memberfunction:`(7)`: returns a reference to the basis object

   - Note: do NOT use ``const_cast`` on the returned value to modify
     the referenced object, as doing so will break the invariant
     and the immutability of the class.
     If you need a different linear subspace, create a new object.

:memberfunction:`(8)`: returns a reference to the offset object

   - Same "rule" about avoiding const casting applies here as in (7) above

:memberfunction:`(9)`: query if the subspace is affine or not


Notes
^^^^^

- the design and semantics of the constructors and special member
  functions needs a clarification: the intent is to ensure
  the class *always* has value semantics because it nees to work
  even for template arguments with view semantics.
  Therefore, anywhere we have an lvalue argument, we cannot
  simply use a copy-constructor to construct the corresponding
  data member because that would break the semantics.
  Suppose that ``basis`` were an instance of a ``Tpetra::MultiVector``
  or a ``Kokkos::View``, which are both classes with view semantics.
  In such case, if constructor (1) were to simply copy-construct
  ``basis``, it would lead to a shallow copy. This would
  compromise the target immutability and invariant of the class,
  since any changes to the original object would be reflected
  in the data member copy-constructed from it.
  Therefore, by relying on ``pressio::ops::clone()``, we can
  ensure the correct semantics and achieve the desired immutability.

Example
^^^^^^^

.. code-block:: cpp

   #include "pressio/rom_subspaces.hpp"
   int main(){
     using basis_t  = Eigen::MatrixXd;
     using offset_t = Eigen::VectorXd;
     using space_t  = pressio::rom::AffineSubspace<basis_t, offset_t>;

     basis_t columnBasis1(56,5);
     offset_t offset1(56);
     // fill columnBasis1 and offset1
     // construct space via constructor (1)
     space_t space(columnBasis1, offset1, true);
     const auto & basis = space.viewBasis();
     // basis is a reference to an object that is NOT
     // the same as columnBasis1 but contains the same values

     basis_t columnBasis2(56,5);
     // fill columnBasis2
     // construct space via constructor (2)
     space_t space(std::move(columnBasis2), offset_t(56), false);

     basis_t columnBasis3(56,5);
     offset_t offset3(56);
     // fill columnBasis3 and offset3
     // construct space via constructor (3)
     space_t space(columnBasis3, std::move(offset3), false);

   }
