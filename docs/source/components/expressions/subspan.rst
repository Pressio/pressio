.. include:: ../../mydefs.rst

``subspan``
===========

Header: ``<pressio/expressions.hpp>``

API
---

.. code-block:: cpp

  namespace pressio {

  template<class T>
  /*impl defined*/ subspan(T & operand,
			   std::pair<std::size_t, std::size_t> rowsRange,
			   std::pair<std::size_t, std::size_t> colsRange));

  } // end namespace pressio

Parameters
----------

* ``operand``: the object to construct a subspan of

* ``rowsRange``: identifies the rows interval to use (the second index is exclusive)

* ``colsRange``: identifies the cols interval to use (the second index is exclusive)

Constraints
~~~~~~~~~~~

- ``T`` must be:

  - an Eigen dense matrix, ``pressio::is_dense_matrix_eigen<T>::value == true``

  - or a Kokkos rank-2 view, i.e. ``pressio::is_dense_matrix_kokkos<T>::value == true``

Preconditions
~~~~~~~~~~~~~

:red:`finish`

Mandates
~~~~~~~~

:red:`finish`


Return value
~~~~~~~~~~~~

None

Effects
~~~~~~~

Returns an expression object that represents a subspan of the ``operand``.
We refer to this as an expression because the subspan function does not allocate
new memory or copy data, but only creates an instance of a
class that "represents the subspan operation".

Postconditions
~~~~~~~~~~~~~~

The returned object is valid to use until the operand goes out of scope.
If the operand goes out of scope but you still have a subspan object,
the state of the subspan object is undefined.

:red:`finish`
