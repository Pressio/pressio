/*
//@HEADER
// ************************************************************************
//
//                                  at.hpp                                 
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


#ifndef TINYMPL_VARIADIC_AT_HPP
#define TINYMPL_VARIADIC_AT_HPP

#include <cstddef>
#include <cstdlib>
#include <type_traits>
#include <utility>

#include "../sequence.hpp"
#include "../identity.hpp"
#include "fill_n.hpp"

#include <future>
#include <condition_variable>

namespace tinympl {
namespace variadic {

// Constant time at<...> implementation borrowed from the brigand project (search on github)
namespace detail {

template <typename T> struct _element_at;

template <std::size_t... Idxs>
struct _element_at<std::integer_sequence<std::size_t, Idxs...>> {
  template <class T>
  static identity<T> at(
    typename tinympl::ignore_value_argument<std::size_t, Idxs, void const*>::type...,
    identity<T>*,
    ...
  );
};

template <std::size_t i, typename... Args>
struct _at_impl
#if 0
  : decltype(
      // Since std::make_index_sequence is implemented in constant "time" via
      // a compiler extension in most implementations, this implementation of at
      // operates in constant (meta-)time also
      _element_at<std::make_index_sequence<i>>::at(
        static_cast<identity<Args>*>(nullptr)...
      )
    )
#else
  : identity<
      std::tuple_element_t<i, std::tuple<Args...>>
    >
#endif
{ };


} // end namespace detail

/**
 * \ingroup VarBasics
 * \class at
 * \brief Extract the i-th element of a variadic template
 * \param i The index to extract
 */
template <std::size_t i, class... Args>
using at = detail::_at_impl<i, Args...>;

template<std::size_t i, class... Args>
using at_t = typename at<i, Args...>::type;


namespace types_only {

template <typename WrappedSpot, class... Args>
struct at : public tinympl::variadic::at<WrappedSpot::value, Args...>
{ };

} // end namespace types_only

template <typename Default, std::size_t i, typename... Args>
struct at_or : tinympl::identity<
  typename std::conditional_t<
    i < sizeof...(Args),
    tinympl::variadic::at<i, Args...>,
    tinympl::identity<Default>
  >::type
> { };

template <typename Default, std::size_t i, typename... Args>
using at_or_t = typename at_or<Default, i, Args...>::type;


} // namespace variadic
} // namespace tinympl

#endif // TINYMPL_VARIADIC_AT_HPP
