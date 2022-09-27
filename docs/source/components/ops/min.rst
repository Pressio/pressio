.. include:: ../../mydefs.rst

``min``
==========

Header: ``<pressio/ops.hpp>``

API
---

.. code-block:: cpp

  namespace pressio { namespace ops{

  template <typename T>
  min(const T & obj);

  }} // end namespace pressio::ops

Parameters
----------

* ``obj``: the object to be read from

Constraints
~~~~~~~~~~~

- ``T`` must be:

  - an Eigen vector or matrix object, i.e. one of the following is
    true ``pressio::is_vector_eigen<T>::value``, ``pressio::is_dense_matrix_eigen<T>::value``,
    ``pressio::is_sparse_matrix_eigen<T>::value`` or

  - a Pressio expression, e.g., span, diag, subspan, operating on an Eigen object.

Mandates
~~~~~~~~

:red:`finish`

Preconditions
~~~~~~~~~~~~~

:red:`finish`

Return value
~~~~~~~~~~~~

- Returns the minimum of all coefficients.

Effects
~~~~~~~

- Calculates and returns the minimum of all coefficients.

- This is a blocking operation, i.e. the kernel completes before returning.

Postconditions
~~~~~~~~~~~~~~

:red:`finish`

