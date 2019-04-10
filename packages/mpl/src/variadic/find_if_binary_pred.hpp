#ifndef ROMPP_MPL_VARIADIC_FIND_IF_BINARY_HPP
#define ROMPP_MPL_VARIADIC_FIND_IF_BINARY_HPP

#include <type_traits>
#include <cstddef>

namespace rompp{ namespace mpl{ namespace variadic {

/**
 * \ingroup VarNonModAlgs
 * \class find_if_binary_pred Compute the index of the first element in the sequence Args...
 which satisfies a given binary predicate

 * \param Predicate The test predicate - `F<T, attribute_t>::type::value` shall be convertible to bool

 * \param attribute_t the second template parameter for the binary predicate

 * \param Args... the input sequence to test

 * \return `find_if<...>::type` is `std::integral_constant<std::size_t,v>` where
`v` is the 0-based index of the first element which satisfy `F`. If no such
element exists, `v` is `size<Sequence>::value`.

 * \sa tinympl::find_if_binary_pred
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

#endif // ROMPP_MPL_VARIADIC_FIND_IF_HPP
