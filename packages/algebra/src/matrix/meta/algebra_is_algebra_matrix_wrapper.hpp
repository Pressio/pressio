
#ifndef ALGEBRA_IS_ALGEBRA_MATRIX_WRAPPER_HPP_
#define ALGEBRA_IS_ALGEBRA_MATRIX_WRAPPER_HPP_

#include "../algebra_matrix_traits.hpp"

namespace rompp{ namespace algebra{ namespace meta {

template <typename T, typename enable = void>
struct is_algebra_matrix_wrapper : std::false_type {};

template <typename T>
struct is_algebra_matrix_wrapper< T,
		       typename
		       std::enable_if<
			 algebra::details::traits<T>::is_matrix
			 >::type
		       > : std::true_type{};


#define STATIC_ASSERT_IS_ALGEBRA_MATRIX_WRAPPER(TYPE) \
  static_assert( algebra::meta::is_algebra_matrix_wrapper<TYPE>::value, \
		 "THIS_IS_NOT_A_ALGEBRA_MATRIX_WRAPPER")

}}}//end namespace rompp::algebra::meta
#endif
