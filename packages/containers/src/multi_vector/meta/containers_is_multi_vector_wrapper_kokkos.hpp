
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_IS_MULTI_VECTOR_WRAPPER_KOKKOS_HPP_
#define CONTAINERS_IS_MULTI_VECTOR_WRAPPER_KOKKOS_HPP_

#include "../containers_multi_vector_traits.hpp"

namespace pressio{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_multi_vector_wrapper_kokkos : std::false_type {};

template <typename T>
struct is_multi_vector_wrapper_kokkos<
  T, ::pressio::mpl::enable_if_t<
       containers::details::traits<T>::is_multi_vector &&
       containers::details::traits<T>::wrapped_multi_vector_identifier==
       containers::details::WrappedMultiVectorIdentifier::Kokkos
       >
  >
  : std::true_type{};

}}}//end namespace pressio::containers::meta
#endif
#endif
