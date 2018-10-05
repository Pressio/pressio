/*
//@HEADER
// ************************************************************************
//
//                      is_metafunction.hpp
//                          tinympl
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

#ifndef TINYMPL_IS_METAFUNCTION_HPP
#define TINYMPL_IS_METAFUNCTION_HPP

#include "detection.hpp"
#include "bool.hpp"

namespace tinympl {

namespace _impl {

template <typename T>
using _has_nullary_apply_template_archetype = typename T::template apply<>;

template <typename T>
using _has_nullary_apply_template =
  is_detected< _has_nullary_apply_template_archetype, T >;

template <typename T>
using _has_unary_apply_template_archetype = typename T::template apply<void>;

template <typename T>
using _has_unary_apply_template =
  is_detected< _has_unary_apply_template_archetype, T >;

template <typename T>
using _has_binary_apply_template_archetype =
  typename T::template apply<void, void>;

template <typename T>
using _has_binary_apply_template =
  is_detected< _has_binary_apply_template_archetype, T >;

template <typename T>
using _has_64ary_apply_template_archetype =
  // Note that this doesn't isn't useful if someone happens to have a non-variadic
  // template with exactly 64 arguments, but in that case they're probably trying
  // to emulate variadic arguments anyway
  typename T::template apply<
    void, void, void, void, void, void, void, void, void, void,
    void, void, void, void, void, void, void, void, void, void,
    void, void, void, void, void, void, void, void, void, void,
    void, void, void, void, void, void, void, void, void, void,
    void, void, void, void, void, void, void, void, void, void,
    void, void, void, void, void, void, void, void, void, void,
    void, void, void, void
  >;

template <typename T>
using _has_64ary_apply_template =
  is_detected< _has_64ary_apply_template_archetype, T >;

template <typename T, typename Enable=void>
struct _has_apply_template
  : std::false_type
{ };

template <typename T>
struct _has_apply_template<T,
  void_template_t<T::template apply>
> : std::true_type
{ };

} // end namespace _impl

template <typename T>
struct is_variadic_metafunction_class
  : bool_<
    // We can't say for sure, but if it accepts 1, 2, and 64 args, it's probably
    // variadic (or at least a valiant attempt to emulate variadic args)
    _impl::_has_nullary_apply_template<T>::value
      and _impl::_has_unary_apply_template<T>::value
      and _impl::_has_binary_apply_template<T>::value
      and _impl::_has_64ary_apply_template<T>::value
  >
{ };

template <typename T>
struct is_nullary_metafunction_class
  : bool_< _impl::_has_nullary_apply_template<T>::value >
{ };

template <typename T>
struct is_binary_metafunction_class
  : bool_< _impl::_has_binary_apply_template<T>::value >
{ };

template <typename T>
struct is_exclusively_unary_metafunction_class
  : bool_<
      _impl::_has_unary_apply_template<T>::value
      and not _impl::_has_binary_apply_template<T>::value
      and not _impl::_has_nullary_apply_template<T>::value
    >
{ };

template <typename T>
struct is_exclusively_binary_metafunction_class
  : bool_<
      _impl::_has_binary_apply_template<T>::value
      and not _impl::_has_unary_apply_template<T>::value
      and not _impl::_has_64ary_apply_template<T>::value
    >
{ };

template <typename T>
struct is_metafunction_class
  : bool_<
      _impl::_has_apply_template<T>::value
    >
{ };

} // end namespace tinympl

#endif //TINYMPL_IS_METAFUNCTION_HPP
