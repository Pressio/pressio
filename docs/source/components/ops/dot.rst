.. include:: ../../mydefs.rst

``dot``
========

Header: ``<pressio/ops.hpp>``

API
---

.. code-block:: cpp

  namespace pressio { namespace ops{

  template<typename vec_type>
  pressio::Traits<vec_type>::scalar_type dot(const vec_type & a, const vec_type & b);

  template<typename vec_type>
  void dot(const vec_type & a, const vec_type & b, typename ::pressio::Traits<vec_type>::scalar_type & result);

  }} // end namespace pressio::ops

Parameters
----------

* ``a`` and ``b``: the objects to be used for calculation of dot product

* [if present] ``result``: the calculated value of dot product

Constraints
~~~~~~~~~~~

- ``vec_type`` must be:

  - an Eigen vector, i.e. ``pressio::is_vector_eigen<vec_type>::value`` or

  - a pressio expression, e.g., span, operating on an Eigen object or

  - a Epetra vector, i.e. ``pressio::is_vector_epetra<vec_type>::value`` or

  - a Kokkos rank-1, i.e. ``pressio::is_vector_kokkos<vec_type>::value`` or

  - a Tpetra vector, i.e. ``pressio::is_vector_tpetra<vec_type>::value`` or

  - a Tpetra block vector, i.e. ``pressio::is_vector_tpetra_block<vec_type>::value``

Mandates
~~~~~~~~

:red:`finish`

Preconditions
~~~~~~~~~~~~~

:red:`finish`

Return value
~~~~~~~~~~~~

- Depending on the overload it will:

  - return the calculated value or

  - return ``void``

Effects
~~~~~~~

- Calculates the dot product of ``a`` and ``b``. The calculated value is of type ``::pressio::Traits<vec_type>::scalar_type``

- Depending on the overload it will:

  - return the calculated value or

  - update the parameter ``result``

Postconditions
~~~~~~~~~~~~~~

:red:`finish`

