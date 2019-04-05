#ifndef TINYMPL_VARIADIC_FIND_IF_HPP
#define TINYMPL_VARIADIC_FIND_IF_HPP

#include <type_traits>
#include <cstddef>

namespace tinympl { namespace variadic {

/**
 * \ingroup VarNonModAlgs
 * \class find_if Compute the index of the first element in the sequence which satisfies a given predicate
 * \param Predicate The test predicate - `F<T>::type::value` shall be convertible to bool
 * \param Args... the input sequence
 * \return `find_if<...>::type` is `std::integral_constant<std::size_t,v>` where
`v` is the 0-based index of the first element which satisfy `F`. If no such
element exists, `v` is `size<Sequence>::value`.
 * \sa tinympl::find_if
 */
template<template<class ...> class Predicate, class ... Args>
struct find_if;

template<template<class ...> class Predicate>
struct find_if<Predicate>
  : std::integral_constant<std::size_t, 0>
{};

template<template<class ...T> class Predicate, class Head, class ... Tail>
struct find_if<Predicate, Head, Tail...>
  : std::conditional <
  Predicate<Head>::type::value,
  std::integral_constant<std::size_t, 0>,
  std::integral_constant <
    std::size_t, 1 + find_if<Predicate, Tail...>::type::value
    >
  >::type
{};

template <template <class... T> class Predicate, class... Args>
using find_if_t = typename find_if<Predicate, Args...>::type;

template <template <class... T> class Predicate, class... Args>
constexpr auto find_if_v = find_if<Predicate, Args...>::type::value;

}} // namespace variadic::tinympl

#endif // TINYMPL_VARIADIC_FIND_IF_HPP
