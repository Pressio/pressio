.. include:: ../mydefs.rst

.. role:: cpp(code)
   :language: cpp

Linear Subspace
===============

Header: ``<pressio/rom_subspaces.hpp>``

Synopsis
--------

.. code-block:: cpp

   namespace pressio{ namespace rom{

   template <class BasisMatrixType>
   class LinearSubspace
   {
     public:
       enum class SpanningSet{Columns, Rows};
       using basis_matrix_type = /*see description*/;

       // constructors
       LinearSubspace(const basis_matrix_type & basisMatrix,            (1)
                      SpanningSet spanningSetValue);

       LinearSubspace(basis_matrix_type && basisMatrix,                 (2)
                      SpanningSet spanningSetValue);

       // copy constructor and copy assignment
       LinearSubspace(const LinearSubspace &);                          (3)
       LinearSubspace & operator=(const LinearSubspace&) = delete;      (4)

       // accessors
       const basis_matrix_type & basis() const;                         (5)
       std::size_t dimension() const;                                   (6)
       bool isColumnSpace() const;                                      (7)
       bool isRowSpace() const;                                         (8)

     private:
       const basis_matrix_type basis_; // exposition-only
   };

   }} //end namespace

Constraints
-----------

- :cpp:`std::is_copy_constructible<basis_matrix_type>::value == true` \
  :cpp:`&& std::is_pointer<basis_matrix_type>::value == false` \
  :cpp:`&& pressio::mpl::is_std_shared_ptr<basis_matrix_type>::value == false`

Member types
------------

- ``basis_matrix_type``: same as ``std::remove_cv_t<BasisMatrixType>``


Member functions
----------------

Constructors
^^^^^^^^^^^^

:memberfunction:`(1)`: constructor accepting an lvalue

  *Parameters*:

  - ``basisMatrix``: the basis matrix object

  - ``spanningSetValue``: value from ``SpanningSet`` to specify whether
    ``basisMatrix`` is a row or column span

  *Preconditions*:

  - ``basisMatrix`` must *represent a linear basis* such that

    - if ``spanningSetValue == SpanningSet::Columns``, ``basisMatrix`` should be *full column rank*

    - if ``spanningSetValue == SpanningSet::Rows``, ``basisMatrix`` should be *full row rank*

  - the following operation ``auto basisMatrixClone = pressio::ops::clone(basisMatrix)``
    conforms to the requirements/semantics as `explained here <ops/clone.html>`__;
    in brief, ``basisMatrixClone`` is an object independent of ``basisMatrix`` and if any
    operation is applied to ``basisMatrixClone`` (or any data it references/uses), they
    do not affect the original ``basisMatrix`` object (or any data it references/uses),
    and viceversa

  *Effects*:

  - makes a new allocation to construct ``basis_`` by invoking the
    copy-constructor as ``basis_(pressio::ops::clone(basisMatrix))``

  *Postconditions*:

  - the constructed object is immutable and ``basis_`` is an
    identical (but distinct) copy of the object ``basisMatrix`` originally
    passed to the constructor. Therefore, this class
    has a clear invariant that is established upon construction.


:memberfunction:`(2)`: constructor accepting an rvalue

  *Parameters*:

  - ``basisMatrix``: the basis matrix object

  - ``spanningSetValue``: value from ``SpanningSet`` to specify whether
    ``basisMatrix`` is a row or column span

  *Preconditions*:

  - ``basisMatrix`` must *represent a linear basis* such that

    - if ``spanningSetValue == SpanningSet::Columns``, ``basisMatrix`` should be *full column rank*

    - if ``spanningSetValue == SpanningSet::Rows``, ``basisMatrix`` should be *full row rank*

  *Effects*:

  - tries to initialize ``basis_`` via move semantics by moving ``basisMatrix``.
    If the type ``basis_matrix_type`` implements move semantics,
    the temporary object is moved-constructed into ``basis_`` and no new
    memory allocation should occur.
    If move semantics are not implemented and the move operation
    falls back to a copy, a new allocation (potentially) occurs.

  *Postconditions*:

  - the constructed object is immutable and ``basis_`` is either
    move constructed from the
    original argument, or is just a copy. As in (1), this class
    has a clear invariant that is established upon construction.

:memberfunction:`(3)`: copy constructor

  *Parameters*:

  - ``other``: another subspace object to be used as source to copy from

  *Effects*:

  - initializes ``basis_`` by calling ``pressio::ops::clone(other.basis_)``,
    therefore the semantics are similar to constructor (1)

  *Postconditions*:

  - the new object will be identical to ``other``


:memberfunction:`(4)`: copy assignment operator

   - this is explicitly marked as deleted because it would violate immutability

Accessors
^^^^^^^^^

:memberfunction:`(5)`: returns a reference to the basis

   - Note: do NOT use ``const_cast`` on the returned value to modify
     the referenced object, as doing so will break the invariant
     and the immutability of the class.
     If you need a different linear subspace, create a new object.

:memberfunction:`(6)`: returns the dimension of the subspace

  *Effects*: let ``S`` be an instance of the class and ``M`` be the basis matrix
  stored in ``S``, then the method returns ``pressio::ops::extent(basisMatrix, i)``
  where ``i==0`` for a row subspace, and ``i==1`` for a column subspace.

:memberfunction:`(7)`: query if the spanning set is the set of columns

  *Effects*: ``true`` for a column subspace, ``false`` otherwise

:memberfunction:`(8)`: query if the spanning set is the set of rows

  *Effects*: ``true`` for a row subspace, ``false`` otherwise


Notes
^^^^^

- the immutability and invariance of this class strictly
  hinges on the assumption that users do not violate const
  correctness: please do not const_cast the return of (5)


- the design and semantics of the class functions needs a clarification:
  the intent is to ensure the class *always* has value semantics because
  it nees to work even for template arguments with view semantics.
  Suppose that ``basis`` were an instance of a ``Tpetra::MultiVector``
  or a ``Kokkos::View``, which are both classes with view semantics.
  In such case, if constructor (1) were to simply copy-construct
  ``basis``, it would lead to a shallow copy. This would
  compromise the immutability and invariant of the class,
  since any changes to the original object would be reflected
  in the data member copy-constructed from it.
  Therefore, in constructor (1) and in (3), we rely
  on ``pressio::ops::clone()``, to ensure the correct semantics
  and achieve the desired immutability.

Example
^^^^^^^

.. code-block:: cpp

   int main(){
     using basis_t  = Eigen::MatrixXd;
     using space_t  = pressio::rom::LinearSubspace<basis_t>;

     basis_t columnBasis1(56,5);
     // fill such that columns are linearly independent
     space_t space(columnBasis1, space_t::SpanningSet::Columns);
     // space.dimension() == 5

     basis_t columnBasis2(56,5);
     // fill such that columns are linearly independent
     space_t space(std::move(columnBasis2), space_t::SpanningSet::Columns);
     // space.dimension() == 5

     basis_t rowBasis1(3,56);
     // fill such that rows are linearly independent
     space_t space(rowBasis1, space_t::SpanningSet::Rows);
     // space.dimension() == 3
   }




..
   - we intentiatlly chose the name ``BasisMatrixType`` rather than limiting it in scope,
     for example narrowing it to ``basisMatrix``, to be more generic for two
     main reasons:

     1. first, the including the word "matrix" in the name is too constraining
	and tied to an implementation detail, since
	a basis could be as well as abstracted in a different way;

     2. second, this flexibility lets users potentially exploit the power of templates.
	In fact, we envision at least two scenarios: the most basic one
	is when ``basis`` is an instance of a ``basis_matrix_type`` class representing
	a "concrete" matrix (e.g. Eigen::Matrix, Tpetra::MultiVector, etc),
	with its columns or rows containing the basis vectors spanning the space.
	An alternative more flexible (but more complicated to handle) scenario is one
	where ``basis`` is an instance of a class that is NOT a "concrete matrix"
	but an abstraction/wrapper of one. This could be an "expert" usage useful in
	some situations.
