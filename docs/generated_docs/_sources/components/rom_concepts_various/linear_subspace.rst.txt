.. include:: ../../mydefs.rst

``VectorSubspace``
==================

Header: ``<pressio/rom_concepts.hpp>``

.. literalinclude:: ../../../../include/pressio/rom/concepts/linear_subspace.hpp
   :language: cpp
   :lines: 52-79

Semantic requirements
---------------------

Given an instance ``A`` of type ``T``, ``VectorSubspace<T>``
is modeled if it is satisfied, all subsumed concepts are modeled and:

- ``A`` represents a subspace closed under addition and scalar multiplication

- upon construction, ``A`` becomes immutable and,
  consequently, the underlying basis and subspace it represents do not change

- the copy constructor ``auto B = A`` must be such that
  ``B`` and ``A`` are independent objects but both
  representing the same subspace. In other words, copying ``A``
  does not have view semantics.

- move semantics should be the same as copy semantics, or it would
  violate the immutability

- let ``auto & basis = A.basis()``, then:

  - if ``A.isColumnSpace() == true``, then ``basis`` is full *column* rank

  - if ``A.isRowSpace() == true``, then ``basis`` is full *row* rank

- if ``A.isColumnSpace() == true``, then it ``A.isRowSpace() == false`` and vice versa

- ``auto dim = A.dimension()`` represents the true dimensionality of the subspace



Syntax-only snippet
-------------------

.. code-block:: cpp

   class SampleClass
   {
     public:
       using basis_matrix_type = Eigen::MatrixXd;

       const basis_matrix_type & basis() const;
       std::size_t dimension() const;
       bool isColumnSpace() const;
       bool isRowSpace() const;

    private:
      /*
        Since we have a const-qualified data member, copy and move assignments
        do not participate in overload resolution, so the class is not
	copy assignable or move assignable which is needed to satisfy the concept.
      */
      const basis_matrix_type matrix_;
   }
