.. include:: ../../mydefs.rst

``abs``
========

Header: ``<pressio/ops.hpp>``

API
---

.. code-block:: cpp

  namespace pressio { namespace ops{

  template <class T1, class T2>
  void abs(T1 & destination, const T2 & source);

  }} // end namespace pressio::ops

Parameters
----------

* ``destination``: the object to write to

* ``source``: the object whose elements are used to compute the absolute value of

Constraints
~~~~~~~~~~~

- ``T1`` and ``T2`` must be:

  - a Pressio expression, e.g., span, operating on an Eigen object or Kokkos rank-1 object or

  - an Eigen vector, i.e. ``pressio::is_vector_eigen<T>::value`` or

  - an Epetra vector, i.e. ``pressio::is_vector_epetra<T>::value`` or

  - a Kokkos rank-1 view, i.e. ``pressio::is_vector_kokkos<T>::value`` or

  - a Tpetra vector, i.e. ``pressio::is_vector_tpetra<T>::value`` or

  - a Tpetra block vector, i.e. ``pressio::is_vector_tpetra_block<T>::value``


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

- Overwrite each element of ``destination`` with abs value of ``source``.

- This is a blocking operation, i.e. the kernel completes before returning.

Postconditions
~~~~~~~~~~~~~~

:red:`finish`
