/*
//@HEADER
// ************************************************************************
//
//                    contains_placeholder.hpp
//                         whatever
//              Copyright (C) 2017 NTESS, LLC
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

#ifndef META_TINYMPL_CONTAINS_PLACEHOLDERS_HPP_
#define META_TINYMPL_CONTAINS_PLACEHOLDERS_HPP_

#include "bind.hpp"
#include "variadic/any_of.hpp"

namespace tinympl {

template <class... T>
struct contains_placeholder
 : std::integral_constant<bool, is_placeholder<T...>::value != 0>
{ };

template <>
struct contains_placeholder<>
 : public std::false_type
{ };


namespace _impl {

// This fixes an error arising from "contains_placeholder" being
// interpreted as a fully-specialized type in the specializations below

template <class T>
struct _contains_placeholder {
  typedef typename contains_placeholder<T>::type type;
  static constexpr const bool value = type::value;
};

}

// an overload for 2 or more variadic arguments

template <class T1, class T2, class... Ts>
struct contains_placeholder<T1, T2, Ts...>
{
  typedef typename variadic::any_of<
    _impl::_contains_placeholder, T1, T2, Ts...
  >::type type;
  static constexpr const bool value = type::value;
};

template <
  template <class...> class F,
  typename... Args
>
struct contains_placeholder<F<Args...>>
{
  typedef typename variadic::any_of<
    _impl::_contains_placeholder, Args...
  >::type type;
  static constexpr const bool value = type::value;
};

} // end namespace tinympl


#endif /* META_TINYMPL_CONTAINS_PLACEHOLDERS_HPP_ */
