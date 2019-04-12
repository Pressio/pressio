
#ifndef ROMPP_MPL_HAS_SUBTRACTION_OPERATOR_HPP_
#define ROMPP_MPL_HAS_SUBTRACTION_OPERATOR_HPP_

#include <type_traits>

namespace rompp{ namespace mpl{

template<typename T, typename U, typename enable = void>
struct has_subtract_op : std::false_type { };

template<typename T, typename U>
struct has_subtract_op<
  T,U,
  typename
  std::enable_if<
    !std::is_void<decltype(std::declval<T>() -
			   std::declval<U>())
		  >::value
    >::type
  > : std::true_type{};

}} // namespace
#endif
