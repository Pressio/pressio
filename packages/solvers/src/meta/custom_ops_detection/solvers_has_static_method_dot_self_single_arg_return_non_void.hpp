/*
//@HEADER
// ************************************************************************
//
// solvers_has_static_method_dot_self_single_arg_return_non_void.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef SOLVERS_SRC_META_SOLVERS_HAS_STATIC_METHOD_DOT_SELF_SINGLE_ARG_RETURN_NON_VOID_HPP_
#define SOLVERS_SRC_META_SOLVERS_HAS_STATIC_METHOD_DOT_SELF_SINGLE_ARG_RETURN_NON_VOID_HPP_

#include <utility>

namespace pressio{ namespace solvers{ namespace meta {

template <
  typename T, typename arg_t, typename result_t,
  typename enable = void>
struct has_static_method_dot_self_single_arg_return_non_void
  : std::false_type{};

template <typename T, typename arg_t, typename result_t>
struct has_static_method_dot_self_single_arg_return_non_void<
  T, arg_t, result_t,
  mpl::enable_if_t<
    std::is_same<
      decltype
      (
       T::template dot_self<result_t>
       (
	std::declval< arg_t const & >()
	)
       ),
      result_t
      >::value
    >
  > : std::true_type{};

}}} //pressio::solvers::meta
#endif