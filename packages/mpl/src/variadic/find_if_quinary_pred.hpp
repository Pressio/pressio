
#ifndef PRESSIO_MPL_VARIADIC_FIND_IF_QUINARY_HPP
#define PRESSIO_MPL_VARIADIC_FIND_IF_QUINARY_HPP

#include <type_traits>
#include <cstddef>

namespace pressio{ namespace mpl{ namespace variadic {

template<typename T1, typename T2, typename T3, typename T4,
	 template<class ...> class Predicate,
	 class ... Args2>
struct find_if_quinary_pred;

template<typename T1, typename T2, typename T3, typename T4,
	 template<class ...> class Predicate>
struct find_if_quinary_pred<T1, T2, T3, T4, Predicate>
  : std::integral_constant<std::size_t, 0>
{};

template<typename T1, typename T2, typename T3, typename T4,
	 template<class ...T> class Predicate,
	 class Head, class ... Tail>
struct find_if_quinary_pred<T1, T2, T3, T4, Predicate, Head, Tail...>
  : std::conditional <
  Predicate<Head, T1,T2,T3,T4>::type::value,
  std::integral_constant<std::size_t, 0>,
  std::integral_constant <
    std::size_t, 1 + find_if_quinary_pred<T1,T2,T3,T4, Predicate, Tail...>::type::value
    >
  >::type
{};

template <typename T1, typename T2, typename T3, typename T4,
	  template <class... T> class Predicate,
	  class... Args>
using find_if_quinary_pred_t = typename find_if_quinary_pred<T1, T2, T3, T4,
							     Predicate,
							     Args...>::type;

}}} // namespace

#endif 
