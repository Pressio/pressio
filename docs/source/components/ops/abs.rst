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

* ``destination``: the object to store absolute value

* ``source``: the object to calculate abs

Constraints
~~~~~~~~~~~

- ``T1`` and ``T2`` must be:

  - an Eigen vector, i.e. ``pressio::is_vector_eigen<T>::value`` or

  - a epetra vector, i.e. ``pressio::is_vector_epetra<T>::value`` or

  - a Kokkos rank-1 view, i.e. ``pressio::is_vector_kokkos<T>::value`` or

  - a tpetra vector, i.e. ``pressio::is_vector_tpetra<T>::value`` or

  - a tpetra block vector, i.e. ``pressio::is_vector_tpetra_block<T>::value``


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

:red:`finish`

Postconditions
~~~~~~~~~~~~~~

:red:`finish`
