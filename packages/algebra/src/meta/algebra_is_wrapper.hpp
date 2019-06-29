
#ifndef ALGEBRA_IS_ALGEBRA_WRAPPER_HPP_
#define ALGEBRA_IS_ALGEBRA_WRAPPER_HPP_

#include "../vector/meta/algebra_is_algebra_vector_wrapper.hpp"
#include "../matrix/meta/algebra_is_algebra_matrix_wrapper.hpp"
#include "../multi_vector/meta/algebra_is_multi_vector_wrapper.hpp"

namespace rompp{ namespace algebra{ namespace meta {

template <typename T, typename enable = void>
struct is_wrapper : std::false_type {};

template <typename T>
struct is_wrapper<
  T,
  typename
  std::enable_if<
    ::rompp::algebra::meta::is_algebra_vector_wrapper<T>::value or
    ::rompp::algebra::meta::is_algebra_multi_vector_wrapper<T>::value or
    ::rompp::algebra::meta::is_algebra_matrix_wrapper<T>::value
    >::type
  > : std::true_type{};

}}}//end namespace rompp::algebra::meta
#endif
