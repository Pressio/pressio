
/*
//@HEADER
// ************************************************************************
//
//                              to_string.hpp                              
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


#ifndef TINYMPL_TO_STRING_HPP
#define TINYMPL_TO_STRING_HPP

#include "string.hpp"

namespace tinympl
{
namespace detail
{

//Handle generic numbers
template<class T,T value,T base = 10,class = void>
struct to_string_impl
{
  typedef typename to_string_impl<T,value / base,base>::type head;
  typedef typename to_string_impl<T,value % base,base>::type tail;
  typedef typename head::template append<tail>::type type;
};

//Handle negative numbers
template<class T, T value, T base>
struct to_string_impl<T,value,base,
  typename std::enable_if<(value < 0)>::type
> {
  typedef typename to_string_impl<T,-value,base>::type tail;
  typedef typename tail::template insert_c<0,'-'>::type type;
};

//Handle one digit numbers
template<class T,T value,T base>
struct to_string_impl<T,value,base,
  typename std::enable_if<(value >= 0 && value < base)>::type
> {
  static_assert( value >= 0 && value < 16,"Base > 16 not supported");

  typedef basic_string<char,
    (value < 10 ?
      '0' + value :
      'a' + value - 10)
  > type;
};

}

/**
 * \ingroup String
 * @{
 */

//! Construct a string from a given integral value of type `T`.
template<class T,T value> using to_string = detail::to_string_impl<T,value>;
template<class T,T value> using to_string_t = typename to_string<T,value>::type;

//! Construct a string from the integer `value`
template<int value> using to_string_i = detail::to_string_impl<int,value>;
template<int value> using to_string_i_t = typename to_string_i<value>::type;

//! Construct a string from the long integer `value`
template<long value> using to_string_l = detail::to_string_impl<long,value>;
template<long value> using to_string_l_t = typename to_string_l<value>::type;

//! Construct a string from the unsigned integer `value`
template<unsigned value> using to_string_u = detail::to_string_impl<unsigned,value>;
template<unsigned value> using to_string_u_t = typename to_string_u<value>::type;

//! Construct a string from the long long integer `value`
template<long long value> using to_string_ll = detail::to_string_impl<long long,value>;
template<long long value> using to_string_ll_t = typename to_string_ll<value>::type;

/** @} */

}

#endif // TINYMPL_TO_STRING_HPP
