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

  - a pressio expression, e.g., span, diag, subspan, operating on an Eigen object or

  - a Kokkos rank-1 or rank-2 view, i.e. ``pressio::is_vector_kokkos<T>::value``, ``pressio::is_dense_matrix_kokkos<T>::value`` or

  - a tpetra vector or multi-vector, i.e. ``pressio::is_vector_tpetra<T>::value``, ``pressio::is_multi_vector_tpetra<T>::value`` or

  - a tpetra block vector or multi-vector, i.e. ``pressio::is_vector_tpetra_block<T>::value``, ``pressio::is_multi_vector_tpetra_block<T>::value``

- ``ScalarType`` must be convertible to ``pressio::Traits<T>::scalar_type``

Mandates
~~~~~~~~

:red:`finish`

Preconditions
~~~~~~~~~~~~~

:red:`finish`

Return value
~~~~~~~~~~~~

None

Effects
~~~~~~~

Overwrite each element of ``operand`` with ``value``.

Postconditions
~~~~~~~~~~~~~~

:red:`finish`
