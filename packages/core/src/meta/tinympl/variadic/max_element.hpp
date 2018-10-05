/*
//@HEADER
// ************************************************************************
//
//                             max_element.hpp                             
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


#ifndef TINYMPL_VARIADIC_MAX_ELEMENT_HPP
#define TINYMPL_VARIADIC_MAX_ELEMENT_HPP

#include "min_element.hpp"
#include "../bind.hpp"

namespace tinympl {
namespace variadic {

/**
 * \ingroup VarMaxMin
 * \class max_element
 * \brief Compute the index of the largest element in a sequence
 * \param Cmp the comparator function; `Cmp<A,B>::type::value` must be
convertible to bool. Defaults to \ref tinympl::less
 * \param Args... the input sequence
 * \return `max_element<...>::type` is an
`std::integral_constant<std::size_t,v>` where `v` is the 0-based index of the
maximum element
 * \sa tinympl::max_element
 */
template<template<class ... > class Cmp, class ... Args>
struct max_element :
    min_element< bind<Cmp, arg2, arg1 >::template eval_t, Args... > {};

} // namespace variadic
} // namespace tinympl

#endif // TINYMPL_VARIADIC_MAX_ELEMENT_HPP
