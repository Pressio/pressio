
#ifndef CORE_NATIVE_STDLIB_VECTOR_META_HPP_
#define CORE_NATIVE_STDLIB_VECTOR_META_HPP_

#include "../core_meta_basic.hpp"
#include <vector>

namespace rompp{ namespace core{ namespace meta {

template <typename T, typename enable = void>
struct is_vector_stdlib : std::false_type {};

template <typename T>
struct is_vector_stdlib<T,
      typename
      std::enable_if<
	std::is_same<T,
	  std::vector<typename T::value_type>
	  >::value and
	// we do not want to have Vector<Vector<...>>
	// so we need to check that the T::value_type is a
	// scalar type or integral type or complex
	(std::is_floating_point<typename T::value_type>::value ||
	 std::is_integral<typename T::value_type>::value ||
	 is_std_complex<typename T::value_type>::value
	 )
	>::type
      > : std::true_type{};

}}}//end namespace rompp::core::meta
#endif
