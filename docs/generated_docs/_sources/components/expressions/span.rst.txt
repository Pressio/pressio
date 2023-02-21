.. include:: ../../mydefs.rst

``span``
========

Header: ``<pressio/expressions.hpp>``

API
---

.. code-block:: cpp

  namespace pressio {

  template<class T>
  /*impl defined*/ span(T & operand,
                        std::pair<std::size_t, std::size_t> indexRange));

  } // end namespace pressio

Parameters
----------

* ``operand``: the object to construct a span of

* ``indexRange``: a std::pair identifying an interval ``[a, b)`` where the second index is exclusive

Constraints
~~~~~~~~~~~

- ``T`` must be:

  - an Eigen vector object: ``pressio::is_vector_eigen<T>::value == true``

  - or a Kokkos rank-1 view, i.e. ``pressio::is_vector_kokkos<T>::value == true``


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

Returns an expression object that represents a target span of the ``operand``.
We refer to this as an expression because the span function does not allocate
new memory but only creates an instance of a class that "represents the span operation".

Postconditions
~~~~~~~~~~~~~~

The returned span object remains valid to use until the operand goes out of scope.
If the operand goes out of scope but you still have a span object, the state of the span
object is undefined.

:red:`finish`
