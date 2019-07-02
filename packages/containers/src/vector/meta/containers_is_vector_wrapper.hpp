
#ifndef CONTAINERS_IS_CONTAINERS_VECTOR_WRAPPER_HPP_
#define CONTAINERS_IS_CONTAINERS_VECTOR_WRAPPER_HPP_

#include "../containers_vector_traits.hpp"

namespace pressio{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_vector_wrapper : std::false_type {};

template <typename T>
struct is_vector_wrapper<T,
	   typename
	   std::enable_if<
	      containers::details::traits<T>::is_vector
	     >::type
	   > : std::true_type{};

#define STATIC_ASSERT_IS_CONTAINERS_VECTOR_WRAPPER(TYPE) \
  static_assert( containers::meta::is_vector_wrapper<TYPE>::value, \
		 "THIS_IS_NOT_A_CONTAINERS_VECTOR_WRAPPER")

}}}//end namespace pressio::containers::meta
#endif
