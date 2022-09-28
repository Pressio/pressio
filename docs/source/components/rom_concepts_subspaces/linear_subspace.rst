.. include:: ../../mydefs.rst

``VectorSubspace``
==================

Header: ``<pressio/rom_concepts.hpp>``

.. literalinclude:: ../../../../include/pressio/rom/concepts/linear_subspace.hpp
   :language: cpp
   :lines: 54-69

Semantic requirements
---------------------

:red:`finish`

..
   The concept is modeled only if it is satisfied,
   all subsumed concepts are modeled and given
   an instance ``A``, of type ``T``, all of the following hold:

   - immutability: ``A`` becomes immutable upon construction
     and, consequently, so are the underlying basis
     and subspace it represents

   - the copy constructor ``auto B = A`` must be such that
     ``B`` and ``A`` are independent objects but both
     representing the same subspace

   - move semantics should be the same copy semantics

   - let ``auto & basis = A.basis()``, then:

     - if ``A.isColumnSpace() == true``, then ``basis`` is full *column* rank

     - if ``A.isRowSpace() == true``, then ``basis`` is full *row* rank

   - if ``A.isColumnSpace() == true``, then it ``A.isRowSpace() == false`` and vice versa

   - ``auto dim = A.dimension()`` represents the true dimensionality of the subspace

   - must be closed under addition and scalar multiplication

   Syntax-only snippet
   -------------------

   The following shows an example class that *satisfies* the concept:

   .. code-block:: cpp

      class SampleClass
      {
	public:
	  /*By having a user-declared copy assignment, the compiler does
	    not generate automatically a move constructor/move assignment,
	    which means that only the copy constr/copy assignment participate
	    in overload resolution, which means we can achieve what we want
	    by simply not declaring move cnstr/assign.*/

	  SampleClass(const SampleClass & other) = default;
	  SampleClass& operator=(const SampleClass & /*other*/) = delete;

	  using basis_matrix_type = Eigen::MatrixXd;

	  const basis_matrix_type & basis() const;
	  std::size_t dimension() const;
	  bool isColumnSpace() const;
	  bool isRowSpace() const;
      }
