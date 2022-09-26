.. include:: ../../mydefs.rst

``add_to_diagonal``
===================

Header: ``<pressio/ops.hpp>``

API
---

.. code-block:: cpp

  namespace pressio { namespace ops{

  template <class T, class ScalarType>
  void add_to_diagonal(T & out, const ScalarType & value);

  }} // end namespace pressio::ops

Parameters
----------

* ``out``: the object to write to

* ``value``: the value to be added to diagonal values of ``out``

Constraints
~~~~~~~~~~~

- ``T`` must be:

  - an Eigen matrix object, ``pressio::is_dense_matrix_eigen<T>::value`` or
    ``pressio::is_sparse_matrix_eigen<T>::value``

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

- Adds value of ``value`` parameter to each diagonal field in the matrix.

- This is a blocking operation, i.e. the kernel completes before returning.

Postconditions
~~~~~~~~~~~~~~

:red:`finish`
