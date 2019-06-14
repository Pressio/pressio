
#ifndef CORE_IS_CORE_WRAPPER_HPP_
#define CORE_IS_CORE_WRAPPER_HPP_

#include "./vector/meta/core_is_core_vector_wrapper.hpp"
#include "./matrix/meta/core_is_core_matrix_wrapper.hpp"
#include "./multi_vector/meta/core_is_multi_vector_wrapper.hpp"

namespace rompp{ namespace core{ namespace meta {

template <typename T, typename enable = void>
struct is_core_wrapper : std::false_type {};

template <typename T>
struct is_core_wrapper<
  T,
  typename
  std::enable_if<
    ::rompp::core::meta::is_core_vector_wrapper<T>::value or
    ::rompp::core::meta::is_core_multi_vector_wrapper<T>::value or
    ::rompp::core::meta::is_core_matrix_wrapper<T>::value
    >::type
  > : std::true_type{};

}}}//end namespace rompp::core::meta
#endif
