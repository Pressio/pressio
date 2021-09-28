
# type_traits


@m_class{m-note m-default}

@parblock
Defined in header: `<pressio/type_traits.hpp>`

Public namespace: `pressio`
@endparblock

## Overview

Provides functionalities for type support and detection.
One of the main design features of pressio is that it supports arbitrary
types via generic programming and type introspection, but also
provides special support for some data types commonly used.

<br/>

## Traits class

\todo: finish

One of the most important things inside `type_traits` is the `Traits` class:

```cpp
namespace pressio{
  template<class T, class = void> struct Traits;
}
```

To understand the purpose and usage of the traits pattern in C++ there are several resources online.
Quoting Bjarne Stroustrup: "Think of a trait as a small object whose main purpose
is to carry information used by another object or algorithm
to determine "policy" or "implementation details".
Pressio uses specializations of this class to gather *in a uniform way*
compile-time information enabling it to reason about types.
The key point here is that *different TPLs use a variety of naming conventions
for nested typedefs and related things*, so there is not easy way to access
similar information from types of various libraries.
This is what motivated us to implement this `type_traits` component.
We need a standard, uniform way to query types for compile-time information.
We currently have traits specialized for types of a few TPLs, like Trilinos, Kokkos, Eigen.
An example of [one such specialization](https://github.com/Pressio/pressio/blob/main/include/pressio/type_traits/traits_vector.hpp) (in this case for Eigen) is:

```cpp
template <typename T>
struct Traits<
  T, ::pressio::mpl::enable_if_t<is_dynamic_vector_eigen<T>::value>
  >
{

  static constexpr bool is_static = false;
  static constexpr bool is_dynamic  = !is_static;
  using scalar_type   = typename T::Scalar;
  using size_type     = typename T::StorageIndex;
  using reference_type = scalar_type &;
  using const_reference_type = scalar_type const &;

  // some other things
};
```

This `Traits` class play a key role when users want to use arbitrary types (i.e. types
which are not known to presso) and to do so, users shoud specialize this class and make
these specialization visibile to pressio to provide information about their generic types. <br/>

For practical examples of how this class is used, see:
- [Newton-Raphson solver](md_pages_components_nonlinsolvers_nr.html)
- [ode explicit steppers](md_pages_components_ode_steppers_explicit.html)
- [ode implicit steppers](md_pages_components_ode_steppers_implicit.html)

<br/>

## Type detection and identification

We support several metafunctions for detecting
data types commonly used from existing TPLs.
The following list is partial, and more will be added as we continue the development.

| Name                                                               | Description                                                                                                                                           |
|--------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------|
| `template<class T>` <br/> `struct is_static_vector_eigen;`         | Static member constant `value` is true <br/> if `T` is a static Eigen vector. <br/> Requires: `PRESSIO_ENABLE_TPL_EIGEN==On`                          |
| `template<class T>` <br/> `struct is_dynamic_vector_eigen;`        | Static member constant `value` is true <br/> if `T` is a dynamic Eigen vector. <br/> Requires: `PRESSIO_ENABLE_TPL_EIGEN==On`                         |
| `template<class T>` <br/> `struct is_vector_eigen;`                | Static member constant `value` is true <br/> if `T` is a static or dynamic Eigen vector. <br/> Requires: `PRESSIO_ENABLE_TPL_EIGEN==On`               |
| `template<class T>` <br/> `struct is_sparse_matrix_eigen;`         | Static member constant `value` is true <br/> if `T` is a static or dynamic sparse Eigen matrix. <br/> Requires: `PRESSIO_ENABLE_TPL_EIGEN==On`        |
| `template<class T>` <br/> `struct is_static_dense_matrix_eigen;`   | Static member constant `value` is true <br/> if `T` is a static dense Eigen matrix. <br/> Requires: `PRESSIO_ENABLE_TPL_EIGEN==On`                    |
| `template<class T>` <br/> `struct is_dynamic_dense_matrix_eigen;`  | Static member constant `value` is true <br/> if `T` is a dynamic dense Eigen matrix. <br/> Requires: `PRESSIO_ENABLE_TPL_EIGEN==On`                   |
| `template<class T>` <br/> `struct is_dense_matrix_eigen;`          | Static member constant `value` is true <br/> if `T` is a static or dynamic dense Eigen matrix. <br/> Requires: `PRESSIO_ENABLE_TPL_EIGEN==On`         |
| `template<class T>` <br/> `struct is_dense_vector_teuchos;`        | Static member constant `value` is true <br/> if `T` is a dense Teuchos vector. <br/> Requires: `PRESSIO_ENABLE_TPL_TRILINOS==On`                      |
| `template<class T>` <br/> `struct is_dense_matrix_teuchos;`        | Static member constant `value` is true <br/> if `T` is a dense Teuchos matrix. <br/> Requires: `PRESSIO_ENABLE_TPL_TRILINOS==On`                      |
| `template<class T>` <br/> `struct is_vector_epetra;`               | Static member constant `value` is true <br/> if `T` is Epetra vector. <br/> Requires: `PRESSIO_ENABLE_TPL_TRILINOS==On`                      |
| `template<class T>` <br/> `struct is_multi_vector_epetra;`         | Static member constant `value` is true <br/> if `T` is an Epetra multi vector. <br/> Requires: `PRESSIO_ENABLE_TPL_TRILINOS==On`                |
| `template<class T>` <br/> `struct is_vector_tpetra;`               | Static member constant `value` is true <br/> if `T` is a Tpetra vector. <br/> Requires: `PRESSIO_ENABLE_TPL_TRILINOS==On`                       |
| `template<class T>` <br/> `struct is_multi_vector_tpetra;`         | Static member constant `value` is true <br/> if `T` is a Tpetra multi vector. <br/> Requires: `PRESSIO_ENABLE_TPL_TRILINOS==On`                 |
| `template<class T>` <br/> `struct is_vector_tpetra_block;`         | Static member constant `value` is true <br/> if `T` is a Tpetra-block vector. <br/> Requires: `PRESSIO_ENABLE_TPL_TRILINOS==On`                 |
| `template<class T>` <br/> `struct is_multi_vector_tpetra_block;`   | Static member constant `value` is true <br/> if `T` is a Tpetra-block multi vector. <br/> Requires: `PRESSIO_ENABLE_TPL_TRILINOS==On`           |
| `template<class T>` <br/> `struct is_static_vector_kokkos;`        | Static member constant `value` is true <br/> if `T` is a static Kokkos vector (rank-1 View). <br/> Requires: `PRESSIO_ENABLE_TPL_KOKKOS==On`                  |
| `template<class T>` <br/> `struct is_dynamic_vector_kokkos;`       | Static member constant `value` is true <br/> if `T` is a dynamic Kokkos vector (rank-1 View).  <br/> Requires: `PRESSIO_ENABLE_TPL_KOKKOS==On`                |
| `template<class T>` <br/> `struct is_vector_kokkos;`               | Static member constant `value` is true <br/> if `T` is a static or dynamic Kokkos vector (rank-1 View).  <br/> Requires: `PRESSIO_ENABLE_TPL_KOKKOS==On`          |
| `template<class T>` <br/> `struct is_static_dense_matrix_kokkos;`  | Static member constant `value` is true <br/> if `T` is a static dense Kokkos matrix (rank-2 View). <br/> Requires: `PRESSIO_ENABLE_TPL_KOKKOS==On`            |
| `template<class T>` <br/> `struct is_dynamic_dense_matrix_kokkos;` | Static member constant `value` is true <br/> if `T` is a dynamic dense Kokkos matrix (rank-2 View). <br/> Requires: `PRESSIO_ENABLE_TPL_KOKKOS==On`           |
| `template<class T>` <br/> `struct is_dense_matrix_kokkos;`         | Static member constant `value` is true <br/> if `T` is a static or dynamic dense Kokkos matrix (rank-2 View).       <br/> Requires: `PRESSIO_ENABLE_TPL_KOKKOS==On` |
