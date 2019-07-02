
#ifndef PRESSIO_MPL_VARIADIC_FIND_IF_BINARY_HPP
#define PRESSIO_MPL_VARIADIC_FIND_IF_BINARY_HPP

#include <type_traits>
#include <cstddef>

namespace pressio{ namespace mpl{ namespace variadic {

/**
 * \param Predicate The test predicate - `F<T, attribute_t>::type::value`
 */
template<typename attribute_t,
	 template<class ...> class Predicate,
	 class ... Args2>
struct find_if_binary_pred;

template<typename attribute_t,
	 template<class ...> class Predicate>
struct find_if_binary_pred<attribute_t, Predicate>
  : std::integral_constant<std::size_t, 0>
{};

template<typename attribute_t,
	 template<class ...T> class Predicate,
	 class Head, class ... Tail>
struct find_if_binary_pred<attribute_t, Predicate, Head, Tail...>
  : std::conditional <
  Predicate<Head, attribute_t>::type::value,
  std::integral_constant<std::size_t, 0>,
  std::integral_constant <
    std::size_t, 1 + find_if_binary_pred<attribute_t, Predicate, Tail...>::type::value
    >
  >::type
{};

template <typename attribute_t,
	  template <class... T> class Predicate,
	  class... Args>
using find_if_binary_pred_t = typename find_if_binary_pred<attribute_t,
						 Predicate,
						 Args...>::type;

}}} // namespace 

#endif
