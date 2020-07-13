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

#ifndef MPL_VARIADIC_SIZE_HPP_
#define MPL_VARIADIC_SIZE_HPP_

namespace pressio{ namespace mpl{ namespace variadic {

/**
 * \class size
 * \brief Compute the size of a variadic template
 * \return `size<Args...>::value` is equivalent to `sizeof ... (Args)`
 */
template<class ... Args> 
struct size;

template<class ... Args> 
using size_t = typename size<Args...>::type;

template<class ... Args> 
struct size : std::integral_constant<std::size_t,sizeof...(Args)>
{};

}}} // namespace pressio::mpl::variadic

#endif  // MPL_VARIADIC_SIZE_HPP_
