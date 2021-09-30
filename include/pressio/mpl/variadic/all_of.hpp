// Copyright (C) 2013, Ennio Barbaro.
//
// Use, modification, and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
// See http://sbabbi.github.io/tinympl for documentation.
//
// You are welcome to contact the author at:
//  enniobarbaro@gmail.com
//

#ifndef MPL_VARIADIC_ALL_OF_HPP_
#define MPL_VARIADIC_ALL_OF_HPP_

namespace pressio{ namespace mpl{ namespace variadic {

/**
 * \ingroup VarNonModAlgs
 * \class all_of
 * \brief Determines whether every element in the sequence satisfies the given
predicate
 * \param F the predicate, `F<T>::value` must be convertible to `bool`
 * \return `all_of<...>::type` is a `std::integral_constant<bool,v>` where `v`
is true iff all the elements in the sequence satisfy the predicate `F`
 * \sa tinympl::all_of
 */
template< template<class ... T> class F, class ... Args>
struct all_of;

template< template<class ... T> class F, class Head, class ... Args>
struct all_of<F, Head, Args...>
  : std::conditional <
  F<Head>::value,
  typename all_of<F, Args...>::type,
  std::integral_constant<bool, false> >::type
{};

template< template<class ... T> class F>
struct all_of<F>
  : std::integral_constant<bool, true>
{};

template< template<class... T> class F, class ... Args>
using all_of_t = typename all_of<F, Args...>::type;

}}}
#endif  // MPL_VARIADIC_ALL_OF_HPP_
