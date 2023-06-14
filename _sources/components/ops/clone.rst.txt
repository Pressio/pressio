.. include:: ../../mydefs.rst

.. role:: cpp(code)
   :language: cpp

``clone``
=========

Header: ``<pressio/ops.hpp>``

API
---

.. code-block:: cpp

  namespace pressio { namespace ops{

  template<class T>
  T clone(const T & operand);

  }} // end namespace pressio::ops

Parameters
----------

* ``operand``: the object to clone

Constraints
~~~~~~~~~~~

- ``T`` must be:

  - an Eigen vector or matrix object:
    ``pressio::is_vector_eigen<T>::value || pressio::is_dense_matrix_eigen<T>::value ||
    pressio::is_sparse_matrix_eigen<T>::value``

  - or a Kokkos rank-1 or rank-2 view: ``pressio::is_vector_kokkos<T>::value ||
    pressio::is_dense_matrix_kokkos<T>::value``

  - or an epetra vector or multi-vector: ``pressio::is_vector_epetra<T>::value ||
    pressio::is_multi_vector_epetra<T>::value``

  - or a tpetra vector or multi-vector: ``pressio::is_vector_tpetra<T>::value ||
    pressio::is_multi_vector_tpetra<T>::value``

  - or a tpetra block vector or multi-vector: ``pressio::is_vector_tpetra_block<T>::value ||
    pressio::is_multi_vector_tpetra_block<T>::value``

Preconditions
~~~~~~~~~~~~~

None

Return value, Effects, and Postconditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Creates and returns a new instance of :cpp:`T` by making a new allocation
  and *copying* all values from :cpp:`operand` into it. So it is an exact
  but independent clone of :cpp:`operand`.

- This is a blocking operation, i.e. the kernel completes before returning.

- This kernel has value semantics, even for types that, by default,
  have view semantics like Kokkos, Tpetra or TpetraBlock.
  This means the following:
  let :cpp:`auto result = clone(operand)`, then any operation applied
  to :cpp:`operand` *after* calling clone will NOT
  have any impact on ``result``.
  And any operation applied to ``result`` will not have any impact on ``operand``.

Notes
~~~~~

- For types that have value semantics, for example ``Eigen::Matrix<...>``,
  Epetra vector or MV, this kernel can be implemented by calling the copy constructor
  and returning the copy

- For Kokkos, Tpetra, or TpetraBlock data types, which by default have view semantics
  (i.e. a copy is a shallow copy), the operation can be implemented by first making a new
  object with extents identical to ``operand``, followed by a deep copy, and then return the result.
