
#ifdef HAVE_KOKKOS
#ifndef CONTAINERS_IS_MATRIX_WRAPPER_KOKKOS_HPP_
#define CONTAINERS_IS_MATRIX_WRAPPER_KOKKOS_HPP_

#include "containers_is_dense_matrix_wrapper_kokkos.hpp"
#include "containers_is_sparse_matrix_wrapper_kokkos.hpp"

namespace pressio{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_matrix_wrapper_kokkos : std::false_type {};

template <typename T>
struct is_matrix_wrapper_kokkos<
  T,
  ::pressio::mpl::enable_if_t<
    is_sparse_matrix_wrapper_kokkos<T>::value or
    is_dense_matrix_wrapper_kokkos<T>::value
    >
  >
  : std::true_type{};

}}}//end namespace pressio::containers::meta
#endif
#endif
