
#ifndef CONTAINERS_IS_CONTAINERS_MATRIX_WRAPPER_HPP_
#define CONTAINERS_IS_CONTAINERS_MATRIX_WRAPPER_HPP_

#include "../containers_matrix_traits.hpp"

namespace pressio{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_matrix_wrapper : std::false_type {};

template <typename T>
struct is_matrix_wrapper< T,
		       typename
		       std::enable_if<
			 containers::details::traits<T>::is_matrix
			 >::type
		       > : std::true_type{};


#define STATIC_ASSERT_IS_CONTAINERS_MATRIX_WRAPPER(TYPE) \
  static_assert( containers::meta::is_matrix_wrapper<TYPE>::value, \
		 "THIS_IS_NOT_A_CONTAINERS_MATRIX_WRAPPER")

}}}//end namespace pressio::containers::meta
#endif
