
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_IS_SPARSE_MATRIX_WRAPPER_KOKKOS_HPP_
#define CONTAINERS_IS_SPARSE_MATRIX_WRAPPER_KOKKOS_HPP_

#include "../containers_matrix_traits.hpp"

namespace pressio{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_sparse_matrix_wrapper_kokkos : std::false_type {};

template <typename T>
struct is_sparse_matrix_wrapper_kokkos<
  T, ::pressio::mpl::enable_if_t<
       containers::details::traits<T>::is_matrix &&
       containers::details::traits<T>::wrapped_matrix_identifier==
       containers::details::WrappedMatrixIdentifier::CrsKokkos
       >
  >
  : std::true_type{};

}}}//end namespace pressio::containers::meta
#endif
#endif
