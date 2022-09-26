.. include:: ../../mydefs.rst

``deep_copy``
=============

Header: ``<pressio/ops.hpp>``

API
---

.. code-block:: cpp

  namespace pressio { namespace ops{

  template<typename T>
  void deep_copy(T & dest, const T & src);

  }} // end namespace pressio::ops

Parameters
----------

* ``dest``: the object to write to

* ``src``: the object to copy the data of

Constraints
~~~~~~~~~~~

- ``T`` must be:

  - an Eigen vector or matrix object, i.e. one of the following is
    true ``pressio::is_vector_eigen<T>::value``, ``pressio::is_dense_matrix_eigen<T>::value``,
    ``pressio::is_sparse_matrix_eigen<T>::value`` or

  - a Pressio expression, e.g., span, diag, subspan, operating on an Eigen object or

  - a Epetra vector or multi-vector, i.e. ``pressio::is_vector_epetra<T>::value``,  ``pressio::is_multi_vector_epetra<T>::value`` or

  - a Kokkos rank-1 or rank-2 view, i.e. ``pressio::is_vector_kokkos<T>::value``, ``pressio::is_dense_matrix_kokkos<T>::value`` or

  - a Tpetra vector or multi-vector, i.e. ``pressio::is_vector_tpetra<T>::value``, ``pressio::is_multi_vector_tpetra<T>::value`` or

  - a Tpetra block vector or multi-vector, i.e. ``pressio::is_vector_tpetra_block<T>::value``, ``pressio::is_multi_vector_tpetra_block<T>::value``

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

- Performs deep copy from ``src`` to  ``dest`` object. Both parameters are not linked between each other.

- This is a blocking operation, i.e. the kernel completes before returning.

Postconditions
~~~~~~~~~~~~~~

:red:`finish`
