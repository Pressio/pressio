.. include:: ../../mydefs.rst

``diag``
========

Header: ``<pressio/expressions.hpp>``

API
---

.. code-block:: cpp

  namespace pressio {

  template<class T>
  /*impl defined*/ diag(T & operand);

  } // end namespace pressio

Parameters
----------

* ``operand``: the object whose diagonal we want to view

Constraints
~~~~~~~~~~~

- ``T`` must be:

  - an Eigen dense matrix, ``pressio::is_dense_matrix_eigen<T>::value == true``

  - or a Kokkos rank-2 view, i.e. ``pressio::is_dense_matrix_kokkos<T>::value == true``

Preconditions
~~~~~~~~~~~~~

- ``operand`` must be a square matrix

Mandates
~~~~~~~~

:red:`finish`


Return value
~~~~~~~~~~~~

None

Effects
~~~~~~~

Returns an expression object that represents the "diagonal" of ``operand``.
We refer to this as an expression because the ``diag`` function does not allocate
new memory or copy any data, but only creates an instance of a
class that "represents the diagonal".

Postconditions
~~~~~~~~~~~~~~

The returned object is valid to use until the operand goes out of scope.

:red:`finish`
