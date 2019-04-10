#ifndef ROMPP_MPL_VARIADIC_FIND_IF_HPP
#define ROMPP_MPL_VARIADIC_FIND_IF_HPP

#include <type_traits>
#include <cstddef>

namespace rompp{ namespace mpl{ namespace variadic {

/**
 * \ingroup VarNonModAlgs

 * \class find_if_unary_pred Compute the index of the first element in the sequence which satisfies a given predicate

 * \param UnaryPredicate The test predicate - `F<T>::type::value` shall be convertible to bool

 * \param Args... the input sequence

 * \return `find_if_unary_pred<...>::type` is `std::integral_constant<std::size_t,v>` where
`v` is the 0-based index of the first element which satisfy `F`. If no such
element exists, `v` is `size<Sequence>::value`.

 * \sa tinympl::find_if_unary_pred
 */

template<template<class ...> class UnaryPredicate, class ... Args>
struct find_if_unary_pred;

template<template<class ...> class UnaryPredicate>
struct find_if_unary_pred<UnaryPredicate>
  : std::integral_constant<std::size_t, 0>
{};

template<template<class ...T> class UnaryPredicate, class Head, class ... Tail>
struct find_if_unary_pred<UnaryPredicate, Head, Tail...>
  : std::conditional <
  UnaryPredicate<Head>::type::value,
  std::integral_constant<std::size_t, 0>,
  std::integral_constant <
    std::size_t, 1 + find_if_unary_pred<UnaryPredicate, Tail...>::type::value
    >
  >::type
{};

template <template <class... T> class UnaryPredicate, class... Args>
using find_if_unary_pred_t = typename find_if_unary_pred<UnaryPredicate, Args...>::type;

}}} // namespace 

#endif // ROMPP_MPL_VARIADIC_FIND_IF_HPP
