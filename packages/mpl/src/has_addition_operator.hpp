
#ifndef PRESSIO_MPL_HAS_ADDITION_OPERATOR_HPP_
#define PRESSIO_MPL_HAS_ADDITION_OPERATOR_HPP_

#include <type_traits>

namespace pressio{ namespace mpl{

template<typename T, typename U=T, typename enable = void>
struct has_addition_op : std::false_type { };

template<typename T, typename U>
struct has_addition_op<
  T,U,
  typename
  std::enable_if<
    !std::is_void<decltype(std::declval<T>() +
			   std::declval<U>())
		  >::value
    >::type
  > : std::true_type{};

}} // namespace
#endif
