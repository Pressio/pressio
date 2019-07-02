
#ifndef PRESSIO_MPL_HAS_SUBSCRIPT_OPERATOR_HPP_
#define PRESSIO_MPL_HAS_SUBSCRIPT_OPERATOR_HPP_

#include <type_traits>

namespace pressio{ namespace mpl{

template<typename T,
	 typename ord_t,
	 typename enable = void>
struct has_subscript_op : std::false_type{};

template<typename T,
	 typename ord_t>
struct has_subscript_op<T,
			ord_t,
			typename
			std::enable_if<
			  !std::is_void<decltype(
						 std::declval<T>()[std::declval<ord_t>()]
						 )
					>::value
			  >::type
			> : std::true_type { };
  // void_t<
  // decltype(std::declval<T>()[std::declval<ord_t>()])
  //   >

}} // namespace
#endif
