
#ifndef CORE_IS_CORE_VECTOR_WRAPPER_HPP_
#define CORE_IS_CORE_VECTOR_WRAPPER_HPP_

#include "../core_vector_traits.hpp"

namespace rompp{ namespace core{ namespace meta {

template <typename T, typename enable = void>
struct is_core_vector_wrapper : std::false_type {};

template <typename T>
struct is_core_vector_wrapper<T,
	   typename
	   std::enable_if<
	      core::details::traits<T>::is_vector
	     >::type
	   > : std::true_type{};

#define STATIC_ASSERT_IS_CORE_VECTOR_WRAPPER(TYPE) \
  static_assert( core::meta::is_core_vector_wrapper<TYPE>::value, \
		 "THIS_IS_NOT_A_CORE_VECTOR_WRAPPER")

}}}//end namespace rompp::core::meta
#endif
