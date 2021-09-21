
# type_traits

Defined in header: `<pressio/type_traits.hpp>`

Public namespace: `pressio`

## Overview

This component includes functionalities for type support and detection.
One of the main design features of pressio is that it supports arbitrary
types via generic programming and type introspection, but also
provides special support for some data types commonly used.

## Type detection and identification

We support several metafunctions for detecting
data types commonly used from existing TPLs.
The following list is partial, and more will be added as we continue the development.

| Name                                                               | Description                                                                                   | Enabled if:                       |
|--------------------------------------------------------------------|-----------------------------------------------------------------------------------------------|-----------------------------------|
| `template<class T>` <br/> `struct is_static_vector_eigen;`         | Self-explanatory                                                                              | `PRESSIO_ENABLE_TPL_EIGEN==On`    |
| `template<class T>` <br/> `struct is_dynamic_vector_eigen;`        | Self-explanatory                                                                              | `PRESSIO_ENABLE_TPL_EIGEN==On`    |
| `template<class T>` <br/> `struct is_vector_eigen;`                | Static member constant `value` is true <br/> if `T` is a static or dynamic vector             | `PRESSIO_ENABLE_TPL_EIGEN==On`    |
| `template<class T>` <br/> `struct is_sparse_matrix_eigen;`         | Self-explanatory                                                                              | `PRESSIO_ENABLE_TPL_EIGEN==On`    |
| `template<class T>` <br/> `struct is_static_dense_matrix_eigen;`   | Self-explanatory                                                                              | `PRESSIO_ENABLE_TPL_EIGEN==On`    |
| `template<class T>` <br/> `struct is_dynamic_dense_matrix_eigen;`  | Self-explanatory                                                                              | `PRESSIO_ENABLE_TPL_EIGEN==On`    |
| `template<class T>` <br/> `struct is_dense_matrix_eigen;`          | Static member constant `value` is true <br/> if `T` is <br/> a static or dynamic dense matrix | `PRESSIO_ENABLE_TPL_EIGEN==On`    |
| `template<class T>` <br/> `struct is_dense_vector_teuchos;`        | Self-explanatory                                                                              | `PRESSIO_ENABLE_TPL_TRILINOS==On` |
| `template<class T>` <br/> `struct is_dense_matrix_teuchos;`        | Self-explanatory                                                                              | `PRESSIO_ENABLE_TPL_TRILINOS==On` |
| `template<class T>` <br/> `struct is_vector_epetra;`               | Self-explanatory                                                                              | `PRESSIO_ENABLE_TPL_TRILINOS==On` |
| `template<class T>` <br/> `struct is_multi_vector_epetra;`         | Self-explanatory                                                                              | `PRESSIO_ENABLE_TPL_TRILINOS==On` |
| `template<class T>` <br/> `struct is_vector_tpetra;`               | Self-explanatory                                                                              | `PRESSIO_ENABLE_TPL_TRILINOS==On` |
| `template<class T>` <br/> `struct is_multi_vector_tpetra;`         | Self-explanatory                                                                              | `PRESSIO_ENABLE_TPL_TRILINOS==On` |
| `template<class T>` <br/> `struct is_vector_tpetra_block;`         | Self-explanatory                                                                              | `PRESSIO_ENABLE_TPL_TRILINOS==On` |
| `template<class T>` <br/> `struct is_multi_vector_tpetra_block;`   | Self-explanatory                                                                              | `PRESSIO_ENABLE_TPL_TRILINOS==On` |
| `template<class T>` <br/> `struct is_static_vector_kokkos;`        | Self-explanatory                                                                              | `PRESSIO_ENABLE_TPL_KOKKOS==On`   |
| `template<class T>` <br/> `struct is_dynamic_vector_kokkos;`       | Self-explanatory                                                                              | `PRESSIO_ENABLE_TPL_KOKKOS==On`   |
| `template<class T>` <br/> `struct is_vector_kokkos;`               | Static member constant `value` is true <br/> if `T` is a static or dynamic vector             | `PRESSIO_ENABLE_TPL_KOKKOS==On`   |
| `template<class T>` <br/> `struct is_static_dense_matrix_kokkos;`  | Self-explanatory                                                                              | `PRESSIO_ENABLE_TPL_KOKKOS==On`   |
| `template<class T>` <br/> `struct is_dynamic_dense_matrix_kokkos;` | Self-explanatory                                                                              | `PRESSIO_ENABLE_TPL_KOKKOS==On`   |
| `template<class T>` <br/> `struct is_dense_matrix_kokkos;`         | Static member constant `value` is true <br/> if `T` is a static or dynamic dense matrix       | `PRESSIO_ENABLE_TPL_KOKKOS==On`   |

## Traits class

```cpp
namespace pressio{
  template<class T, class = void> struct Traits;
}
```

This class contains compile-time information enabling pressio to reason about types.
To understand the purpose and usage of the traits pattern in C++ there are several resources online.
Quoting Bjarne Stroustrup: "Think of a trait as a small object whose main purpose
is to carry information used by another object or algorithm
to determine "policy" or "implementation details".

**need to add much more to this**
