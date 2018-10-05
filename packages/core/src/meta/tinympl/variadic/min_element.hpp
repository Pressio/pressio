/*
//@HEADER
// ************************************************************************
//
//                             min_element.hpp                             
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


#ifndef TINYMPL_VARIADIC_MIN_ELEMENT_HPP
#define TINYMPL_VARIADIC_MIN_ELEMENT_HPP

#include "at.hpp"
#include <type_traits>
#include <cstddef>

namespace tinympl {
namespace variadic {

/**
 * \ingroup VarMaxMin
 * \class min_element
 * \brief Compute the index of the smallest element in a sequence
 * \param Cmp the comparator function; `Cmp<A,B>::type::value` must be
convertible to bool. Defaults to \ref tinympl::less
 * \param Args... the input sequence
 * \return `min_element<...>::type` is an
`std::integral_constant<std::size_t,v>` where `v` is the 0-based index of the
minimum element
 * \sa tinympl::min_element
 */
template<template<class ... > class Cmp, class ... Args> struct min_element;

namespace detail {

template<template<class ...> class Comp, class ... > struct min_element_impl;
template<template<class ...> class Comp, class Head, class ... Tail> struct
min_element_impl<Comp, Head, Tail...> {
  private:
    enum {
      next_min = min_element_impl<Comp, Tail...>::type::value
    };

    enum {
      this_min = ! Comp<at_t<next_min, Tail...>, Head>::type::value
    };

  public:
    typedef std::integral_constant<
      std::size_t,
      ( this_min ? 0 : next_min + 1 )
    > type;
};

template<template<class ... > class Comp, class Head> struct
min_element_impl<Comp, Head> {
  typedef std::integral_constant<std::size_t, 0> type;
};

} // namespace detail

template <template <class...> class Comp, class ... Args>
struct min_element :
  detail::min_element_impl<Comp, Args...>::type
{ };

} // namespace variadic
} // namespace tinympl

#endif // TINYMPL_VARIADIC_MIN_ELEMENT_HPP
