
/*
//@HEADER
// ************************************************************************
//
//                               equal_to.hpp                              
//                         whatever
//              Copyright (C) 2015 Sandia Corporation
// This file was adapted from its original form in the tinympl library.
// The original file bore the following copyright:
//   Copyright (C) 2013, Ennio Barbaro.
// See LEGAL.md for more information.
//
// Under the terms of Contract DE-NA-0003525 with NTESS, LLC,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact somebody@sandia.gov
//
// ************************************************************************
//@HEADER
*/


#ifndef TINYMPL_EQUAL_TO_HPP
#define TINYMPL_EQUAL_TO_HPP

#include <type_traits>

#include "bool.hpp"

namespace tinympl {

/**
 * \ingroup Comparisons
 * \class equal_to
 * \brief Determines whether the types `A` and `B` are equal
 * \return `equal_to<A,B>::type` is a `std::integral_constant<bool,v>` where `v` is true iff `A` and `B` are equal
 * \note The default behaviour is to forward the call to std::is_same. Users are allowed to specialize this metafunction for user-defined types
 */
template<class A,class B> struct equal_to : std::is_same<A,B> {};
template<class T,class U,T t,U u> struct equal_to<
	std::integral_constant<T,t>,
	std::integral_constant<U,u> > : std::integral_constant<bool,t ==u> {};

template <typename ValueType>
struct make_value_equal {
	template <ValueType a, ValueType b>
  struct apply {
		static constexpr auto value = a == b;
		using type = bool_<value>;
	};
	template <ValueType a, ValueType b>
  using apply_t = typename apply<a, b>::type;
};

} // namespace tinympl

#endif // TINYMPL_EQUAL_TO_HPP
