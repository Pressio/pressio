.. include:: ../../mydefs.rst

``fill``
========

Header: ``<pressio/ops.hpp>``

API
---

.. code-block:: cpp

  namespace pressio { namespace ops{

  template<class T, class ScalarType>
  void fill(T & operand, ScalarType const & value);

  }} // end namespace pressio::ops

Parameters
----------

* ``operand``: the object to fill

* ``value``: value to use to fill

Constraints
~~~~~~~~~~~

- ``T`` must be:

  - an Eigen vector or matrix object, i.e. one of the following is
    true ``pressio::is_vector_eigen<T>::value``, ``pressio::is_dense_matrix_eigen<T>::value``,
    ``pressio::is_sparse_matrix_eigen<T>::value`` or

  - a Kokkos rank-1 or rank-2 view, i.e. ``pressio::is_vector_kokkos<T>::value``, ``pressio::is_dense_matrix_kokkos<T>::value`` or

  - a Tpetra vector or multi-vector, i.e. ``pressio::is_vector_tpetra<T>::value``, ``pressio::is_multi_vector_tpetra<T>::value`` or

  - a Tpetra block vector or multi-vector, i.e. ``pressio::is_vector_tpetra_block<T>::value``, ``pressio::is_multi_vector_tpetra_block<T>::value`` or

  - a Epetra vector or multi-vector, i.e. ``pressio::is_vector_epetra<T>::value``, ``pressio::is_multi_vector_epetra<T>::value`` or

  - a pressio expression, i.e. ``pressio::diag``, ``pressio::span``, ``pressio::subspan``, based on Eigen or Kokkkos container

- ``ScalarType`` must be convertible to ``pressio::Traits<T>::scalar_type``

Preconditions
~~~~~~~~~~~~~

None

Return value
~~~~~~~~~~~~

None

Effects
~~~~~~~

Overwrites each element of ``operand`` with ``value``.

Postconditions
~~~~~~~~~~~~~~

Each element of ``operand`` is equal to ``value``.
