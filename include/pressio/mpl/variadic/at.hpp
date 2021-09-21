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

#ifndef MPL_VARIADIC_AT_HPP_
#define MPL_VARIADIC_AT_HPP_

namespace pressio{ namespace mpl{ namespace variadic {

/**
 * \class at
 * \brief Extract the i-th element of a variadic template
 * \param i The index to extract
 */
template <std::size_t i, typename... Args> struct at;

template <std::size_t i, typename... Args>
struct at
  : identity<
  typename std::tuple_element<i, std::tuple<Args...>>::type
  >{};

template<std::size_t i, class... Args>
using at_t = typename at<i, Args...>::type;
//----------------------------------------------------------------

template<typename Default, std::size_t i, typename ... Args>
struct at_or
  : std::conditional<
  i < sizeof ... (Args),
      at<i, Args...>,
      identity<Default>
      >::type
  {};

template<typename Default, std::size_t i, typename ... Args>
using at_or_t = typename at_or<Default, i, Args...>::type;

}}}

#endif  // MPL_VARIADIC_AT_HPP_
