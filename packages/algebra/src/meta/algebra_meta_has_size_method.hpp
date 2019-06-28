
#ifndef ALGEBRA_META_META_HAS_SIZE_METHOD_HPP_
#define ALGEBRA_META_META_HAS_SIZE_METHOD_HPP_

#include <type_traits>

namespace rompp{ namespace algebra{ namespace meta {

template<typename T,
	 typename = void>
struct has_size_method : std::false_type{};

template<typename T>
struct has_size_method<
  T,
  typename
  std::enable_if<
    !std::is_void<
      decltype(std::declval<T>().size())
      >::value
    >::type
  > : std::true_type{};


}}} // namespace rompp::algebra::meta
#endif
