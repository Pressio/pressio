/*
//@HEADER
// ************************************************************************
//
//                               find_if.hpp                               
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


#ifndef TINYMPL_VARIADIC_FIND_IF_HPP
#define TINYMPL_VARIADIC_FIND_IF_HPP

#include <type_traits>
#include <cstddef>

#include "../delay.hpp"
#include "../plus.hpp"

namespace tinympl {
namespace variadic {

/**
 * \ingroup VarNonModAlgs
 * \class find_if Compute the index of a the first element in the sequence which
satisfies a given predicate
 * \param F The test predicate - `F<T>::type::value` shall be convertible to
bool
 * \param Args... the input sequence
 * \return `find_if<...>::type` is `std::integral_constant<std::size_t,v>` where
`v` is the 0-based index of the first element which satisfy `F`. If no such
element exists, `v` is `size<Sequence>::value`.
 * \sa tinympl::find_if
 */
template <template <class... T> class F, class... Args>
struct find_if;

template<template<class ...T> class F, class Head, class ... Tail>
struct find_if<F, Head, Tail...> :
  std::conditional<
    F<Head>::type::value,
    std::integral_constant<std::size_t, 0>,
    delay<
      plus,
      std::integral_constant<std::size_t, 1>,
      find_if<F, Tail...>
    >
  >::type::type
{};

template <template <class... T> class F>
struct find_if<F>
  : std::integral_constant<std::size_t, 0>
{};

template <template <class... T> class F, class... Args>
using find_if_t = typename find_if<F, Args...>::type;

} // namespace variadic
} // namespace tinympl

#endif // TINYMPL_VARIADIC_FIND_IF_HPP
