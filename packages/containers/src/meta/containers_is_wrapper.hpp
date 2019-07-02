
#ifndef CONTAINERS_IS_CONTAINERS_WRAPPER_HPP_
#define CONTAINERS_IS_CONTAINERS_WRAPPER_HPP_

#include "../vector/meta/containers_is_vector_wrapper.hpp"
#include "../matrix/meta/containers_is_matrix_wrapper.hpp"
#include "../multi_vector/meta/containers_is_multi_vector_wrapper.hpp"

namespace pressio{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_wrapper : std::false_type {};

template <typename T>
struct is_wrapper<
  T,
  typename
  std::enable_if<
    ::pressio::containers::meta::is_vector_wrapper<T>::value or
    ::pressio::containers::meta::is_multi_vector_wrapper<T>::value or
    ::pressio::containers::meta::is_matrix_wrapper<T>::value
    >::type
  > : std::true_type{};

}}}//end namespace pressio::containers::meta
#endif
