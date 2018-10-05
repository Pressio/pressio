/*
//@HEADER
// ************************************************************************
//
//                      select.hpp
//                         DARMA
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

#ifndef TINYMPL_SELECT_HPP
#define TINYMPL_SELECT_HPP

#include <type_traits>
#include <tinympl/delay_all.hpp>

namespace tinympl {

template <typename... Args>
struct select_first;

template <typename WrappedBoolean, typename T, typename... Args>
struct select_first<WrappedBoolean, T, Args...>
{
  using type = typename std::conditional_t<
    // Allow delayed evaluation: if there's no value, get ::type::value, else get ::value
    extract_bool_value_potentially_lazy<WrappedBoolean>::value,
    identity<T>,
    select_first<Args...>
  >::type;
};

template <>
struct select_first<>;

template <typename... Args>
using select_first_t = typename select_first<Args...>::type;


namespace _impl {

template <typename T>
struct _select_attorney {
  static constexpr auto this_and_rest_are_false = T::_this_and_rest_are_false;
};

} // end namespace _impl

template <typename... Args>
struct select;

template <typename WrappedBoolean, typename T, typename... Args>
struct select<WrappedBoolean, T, Args...> {
  private:
    // Allow delayed evaluation: if there's no value, get ::type::value, else get ::value
    static constexpr auto _current_value = extract_bool_value_potentially_lazy<WrappedBoolean>::value;

    static constexpr auto _this_and_rest_are_false = not _current_value
      and _impl::_select_attorney<select<Args...>>::this_and_rest_are_false;

    friend struct _impl::_select_attorney<select<WrappedBoolean, T, Args...>>;
  public:


    using type = typename std::conditional_t<
      _current_value, identity<T>, select<Args...>
    >::type;

    static_assert(not _current_value or _impl::_select_attorney<select<Args...>>::this_and_rest_are_false,
      "tinympl::select requires exactly one boolean to evaluate to true"
    );

};

template <>
struct select<> {
  private:
    static constexpr auto _this_and_rest_are_false = true;
    friend struct _impl::_select_attorney<select<>>;
};

template <typename... Args>
using select_t = typename select<Args...>::type;

} // end namespace tinympl


#endif //TINYMPL_SELECT_HPP
