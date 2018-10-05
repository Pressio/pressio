/*
//@HEADER
// ************************************************************************
//
//                                fill_n.hpp                               
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


#ifndef TINYMPL_VARIADIC_FILL_N_HPP
#define TINYMPL_VARIADIC_FILL_N_HPP

#include <cstddef>

namespace tinympl {
namespace variadic {

// O(log(N)) list fill borrowed from brigand:
//   (https://github.com/edouarda/brigand/blob/master/brigand/sequences/filled_list.hpp)
namespace detail {

template<class, class>
struct dup_append_list;

template<template<class...> class List, class... Ts, class... Us>
struct dup_append_list<List<Ts...>, List<Us...>>
{
  using type = List<Ts..., Ts..., Us...>;
};

template<class T, template<class...> class List, std::size_t N>
struct fill_n_impl
  : dup_append_list<
      typename fill_n_impl<T, List, N/2>::type,
      typename fill_n_impl<T, List, N - N/2*2>::type
    >
{ };

template<class T, template<class...> class List>
struct fill_n_impl<T, List, 1> {
  using type = List<T>;
};

template<class T, template<class...> class List>
struct fill_n_impl<T, List, 0> {
  using type = List<>;
};

} // end namespace detail

/**
 * \ingroup VarModAlgs
 * \class fill_n
 * \brief Fills an output sequence with n equal elements
 * \param n The number of elements
 * \param T The type of the elements
 * \param Out The output sequence type
 * \return `fill_n<...>::type` is a type templated from `Out` constructed with n
types equal to `T`
 * \sa tinympl::fill_n
 */
template <std::size_t n, class T, template <class ...> class Out>
using fill_n = detail::fill_n_impl<T, Out, n>;


} // namespace variadic
} // namespace tinympl

#endif // TINYMPL_VARIADIC_FILL_N_HPP
