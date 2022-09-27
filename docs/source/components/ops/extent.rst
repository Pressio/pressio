.. include:: ../../mydefs.rst

``extent``
==========

Header: ``<pressio/ops.hpp>``

API
---

.. code-block:: cpp

  namespace pressio { namespace ops{

  template <typename T, class IndexType>
  ::pressio::Traits<T>::global_ordinal_type extent(const T & oIn, const IndexType i);

  }} // end namespace pressio::ops

Parameters
----------

* ``oIn``: the object to be read from

* ``i``: the dimension index to be read from

Constraints
~~~~~~~~~~~

- ``T`` must be:

  - an Eigen vector or matrix object, i.e. one of the following is
    true ``pressio::is_vector_eigen<T>::value``, ``pressio::is_dense_matrix_eigen<T>::value``,
    ``pressio::is_sparse_matrix_eigen<T>::value`` or

  - a Pressio expression, e.g., span, diag, subspan, operating on an Eigen object or

  - a Kokkos rank-1 or rank-2 view, i.e. ``pressio::is_vector_kokkos<T>::value``, ``pressio::is_dense_matrix_kokkos<T>::value`` or

  - a Epetra vector or multi-vector, i.e. ``pressio::is_vector_epetra<vec_type>::value``,
    ``pressio::is_multi_vector_epetra<T>::value`` or

  - a Tpetra vector or multi-vector, i.e. ``pressio::is_vector_tpetra<T>::value``, ``pressio::is_multi_vector_tpetra<T>::value`` or

  - a Tpetra block vector or multi-vector, i.e. ``pressio::is_vector_tpetra_block<T>::value``, ``pressio::is_multi_vector_tpetra_block<T>::value``

- ``IndexType`` must be convertible to ``pressio::traits::size_type``

Mandates
~~~~~~~~

:red:`finish`

Preconditions
~~~~~~~~~~~~~

:red:`finish`

Return value
~~~~~~~~~~~~

- Returns calculated extent value for passed object.

- The returned value is of type ``pressio::Traits<T>::global_ordinal_type``.

Effects
~~~~~~~

- Calculates and returns extent of passed object.

- This is a blocking operation, i.e. the kernel completes before returning.

Postconditions
~~~~~~~~~~~~~~

:red:`finish`

