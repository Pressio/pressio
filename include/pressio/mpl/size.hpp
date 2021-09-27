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

#ifndef MPL_SIZE_HPP_
#define MPL_SIZE_HPP_

#include "./variadic/size.hpp"

namespace pressio{ namespace mpl{

template<class ... Args>
struct size : variadic::size<Args...>{};

template<class ... Args>
using size_t = typename size<Args...>::type;

}}
#endif  // MPL_VARIADIC_SIZE_HPP_
