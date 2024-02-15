.. role:: raw-html-m2r(raw)
   :format: html

.. include:: ../mydefs.rst

``type_traits``
===============

Header: ``<pressio/type_traits.hpp>``

Public namespace: ``pressio``


Scope
-----

Provides functionalities for type introspection and detection.
One of the main design features of pressio is that it supports arbitrary
types via generic programming and type introspection, but also
provides special support for some data types commonly used.


Traits class
------------

\todo: finish

One of the most important things inside ``type_traits`` is the ``Traits`` class:

.. code-block:: cpp

   namespace pressio{
   template<class T, class = void> struct Traits;
   }

To understand the purpose and usage of the traits pattern in C++ there are several resources online.
Quoting Bjarne Stroustrup: "Think of a trait as a small object whose main purpose
is to carry information used by another object or algorithm
to determine "policy" or "implementation details".
Pressio uses specializations of this class to gather *in a uniform way*
compile-time information enabling it to reason about types.
The key point here is that *different TPLs use a variety of naming conventions
for nested typedefs and related things*\ , so there is not easy way to access
similar information from types of various libraries.
This is what motivated us to implement this ``type_traits`` component.
We need a standard, uniform way to query types for compile-time information.
We currently have traits specialized for types of a few TPLs, like Trilinos, Kokkos, Eigen.
An example of one such specialization (in this case for Eigen) is:

.. code-block:: cpp

   template <typename T>
   struct Traits<
     T, std::enable_if_t<is_dynamic_vector_eigen<T>::value>
     >
   {
     using scalar_type   = typename T::Scalar;
     static constexpr int rank = 1;
   };

This ``Traits`` class play a key role when users want to use arbitrary types (i.e. types
which are not known to presso) and to do so, users shoud specialize this class and make
these specialization visibile to pressio to provide information about their generic types. :raw-html-m2r:`<br/>`

For practical examples of how this class is used, see:

:red:`finish`

..
   * `Newton-Raphson solver <nonlinsolvers_nr.html>`_
   * `ode explicit steppers <ode_steppers_explicit.html>`_
   * `ode implicit steppers <ode_steppers_implicit.html>`_



Type detection and identification
---------------------------------

We support several metafunctions for detecting
data types commonly used from existing TPLs.
The following list is partial, and more will be added as we continue the development.

.. list-table::
   :widths: 45 55
   :header-rows: 1

   * - Name
     - Description
   * - ``template<class T>`` :raw-html-m2r:`<br/>` ``struct is_static_vector_eigen;``
     - Provides static member constant ``value`` equal to ``true`` if ``T`` is a static Eigen vector. :raw-html-m2r:`<br/>` Requires: ``PRESSIO_ENABLE_TPL_EIGEN==On``
   * - ``template<class T>`` :raw-html-m2r:`<br/>` ``struct is_dynamic_vector_eigen;``
     - Provides static member constant ``value`` equal to ``true`` if ``T`` is a dynamic Eigen vector. :raw-html-m2r:`<br/>` Requires: ``PRESSIO_ENABLE_TPL_EIGEN==On``
   * - ``template<class T>`` :raw-html-m2r:`<br/>` ``struct is_vector_eigen;``
     - Provides static member constant ``value`` equal to ``true`` if ``T`` is a static or dynamic Eigen vector. :raw-html-m2r:`<br/>` Requires: ``PRESSIO_ENABLE_TPL_EIGEN==On``
   * - ``template<class T>`` :raw-html-m2r:`<br/>` ``struct is_sparse_matrix_eigen;``
     - Provides static member constant ``value`` equal to ``true`` if ``T`` is a static or dynamic sparse Eigen matrix. :raw-html-m2r:`<br/>` Requires: ``PRESSIO_ENABLE_TPL_EIGEN==On``
   * - ``template<class T>`` :raw-html-m2r:`<br/>` ``struct is_static_dense_matrix_eigen;``
     - Provides static member constant ``value`` equal to ``true`` if ``T`` is a static dense Eigen matrix. :raw-html-m2r:`<br/>` Requires: ``PRESSIO_ENABLE_TPL_EIGEN==On``
   * - ``template<class T>`` :raw-html-m2r:`<br/>` ``struct is_dynamic_dense_matrix_eigen;``
     - Provides static member constant ``value`` equal to ``true`` if ``T`` is a dynamic dense Eigen matrix. :raw-html-m2r:`<br/>` Requires: ``PRESSIO_ENABLE_TPL_EIGEN==On``
   * - ``template<class T>`` :raw-html-m2r:`<br/>` ``struct is_dense_matrix_eigen;``
     - Provides static member constant ``value`` equal to ``true`` if ``T`` is a static or dynamic dense Eigen matrix. :raw-html-m2r:`<br/>` Requires: ``PRESSIO_ENABLE_TPL_EIGEN==On``
   * - ``template<class T>`` :raw-html-m2r:`<br/>` ``struct is_dense_vector_teuchos;``
     - Provides static member constant ``value`` equal to ``true`` if ``T`` is a dense Teuchos vector. :raw-html-m2r:`<br/>` Requires: ``PRESSIO_ENABLE_TPL_TRILINOS==On``
   * - ``template<class T>`` :raw-html-m2r:`<br/>` ``struct is_dense_matrix_teuchos;``
     - Provides static member constant ``value`` equal to ``true`` if ``T`` is a dense Teuchos matrix. :raw-html-m2r:`<br/>` Requires: ``PRESSIO_ENABLE_TPL_TRILINOS==On``
   * - ``template<class T>`` :raw-html-m2r:`<br/>` ``struct is_vector_epetra;``
     - Provides static member constant ``value`` equal to ``true`` if ``T`` is Epetra vector. :raw-html-m2r:`<br/>` Requires: ``PRESSIO_ENABLE_TPL_TRILINOS==On``
   * - ``template<class T>`` :raw-html-m2r:`<br/>` ``struct is_multi_vector_epetra;``
     - Provides static member constant ``value`` equal to ``true`` if ``T`` is an Epetra multi vector. :raw-html-m2r:`<br/>` Requires: ``PRESSIO_ENABLE_TPL_TRILINOS==On``
   * - ``template<class T>`` :raw-html-m2r:`<br/>` ``struct is_vector_tpetra;``
     - Provides static member constant ``value`` equal to ``true`` if ``T`` is a Tpetra vector. :raw-html-m2r:`<br/>` Requires: ``PRESSIO_ENABLE_TPL_TRILINOS==On``
   * - ``template<class T>`` :raw-html-m2r:`<br/>` ``struct is_multi_vector_tpetra;``
     - Provides static member constant ``value`` equal to ``true`` if ``T`` is a Tpetra multi vector. :raw-html-m2r:`<br/>` Requires: ``PRESSIO_ENABLE_TPL_TRILINOS==On``
   * - ``template<class T>`` :raw-html-m2r:`<br/>` ``struct is_vector_tpetra_block;``
     - Provides static member constant ``value`` equal to ``true`` if ``T`` is a Tpetra-block vector. :raw-html-m2r:`<br/>` Requires: ``PRESSIO_ENABLE_TPL_TRILINOS==On``
   * - ``template<class T>`` :raw-html-m2r:`<br/>` ``struct is_multi_vector_tpetra_block;``
     - Provides static member constant ``value`` equal to ``true`` if ``T`` is a Tpetra-block multi vector. :raw-html-m2r:`<br/>` Requires: ``PRESSIO_ENABLE_TPL_TRILINOS==On``
   * - ``template<class T>`` :raw-html-m2r:`<br/>` ``struct is_static_vector_kokkos;``
     - Provides static member constant ``value`` equal to ``true`` if ``T`` is a static Kokkos vector (rank-1 View). :raw-html-m2r:`<br/>` Requires: ``PRESSIO_ENABLE_TPL_KOKKOS==On``
   * - ``template<class T>`` :raw-html-m2r:`<br/>` ``struct is_dynamic_vector_kokkos;``
     - Provides static member constant ``value`` equal to ``true`` if ``T`` is a dynamic Kokkos vector (rank-1 View).  :raw-html-m2r:`<br/>` Requires: ``PRESSIO_ENABLE_TPL_KOKKOS==On``
   * - ``template<class T>`` :raw-html-m2r:`<br/>` ``struct is_vector_kokkos;``
     - Provides static member constant ``value`` equal to ``true`` if ``T`` is a static or dynamic Kokkos vector (rank-1 View).  :raw-html-m2r:`<br/>` Requires: ``PRESSIO_ENABLE_TPL_KOKKOS==On``
   * - ``template<class T>`` :raw-html-m2r:`<br/>` ``struct is_static_dense_matrix_kokkos;``
     - Provides static member constant ``value`` equal to ``true`` if ``T`` is a static dense Kokkos matrix (rank-2 View). :raw-html-m2r:`<br/>` Requires: ``PRESSIO_ENABLE_TPL_KOKKOS==On``
   * - ``template<class T>`` :raw-html-m2r:`<br/>` ``struct is_dynamic_dense_matrix_kokkos;``
     - Provides static member constant ``value`` equal to ``true`` if ``T`` is a dynamic dense Kokkos matrix (rank-2 View). :raw-html-m2r:`<br/>` Requires: ``PRESSIO_ENABLE_TPL_KOKKOS==On``
   * - ``template<class T>`` :raw-html-m2r:`<br/>` ``struct is_dense_matrix_kokkos;``
     - Provides static member constant ``value`` equal to ``true`` if ``T`` is a static or dynamic dense Kokkos matrix (rank-2 View).       :raw-html-m2r:`<br/>` Requires: ``PRESSIO_ENABLE_TPL_KOKKOS==On``
