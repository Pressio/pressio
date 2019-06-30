
#ifndef CONTAINERS_IS_MULTI_VECTOR_WRAPPER_HPP_
#define CONTAINERS_IS_MULTI_VECTOR_WRAPPER_HPP_

#include "../containers_multi_vector_traits.hpp"

namespace rompp{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_multi_vector_wrapper : std::false_type {};

template <typename T>
struct is_multi_vector_wrapper<T,
	   typename
	   std::enable_if<
	     containers::details::traits<T>::is_multi_vector
	     >::type
	   > : std::true_type{};
//------------------------------------------------------------

#define STATIC_ASSERT_IS_CONTAINERS_MULTI_VECTOR_WRAPPER(TYPE) \
  static_assert(containers::meta::is_multi_vector_wrapper<TYPE>::value,\
		"THIS_IS_NOT_A_CONTAINERS_MULTI_VECTOR_WRAPPER")

}}}//end namespace rompp::containers::meta
#endif
