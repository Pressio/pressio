/*
//@HEADER
// ************************************************************************
//
// find_if_quinary_pred.hpp
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

#ifndef MPL_VARIADIC_FIND_IF_QUINARY_PRED_HPP_
#define MPL_VARIADIC_FIND_IF_QUINARY_PRED_HPP_

namespace pressio{ namespace mpl{ namespace variadic {

template<typename T1, typename T2, typename T3, typename T4,
	 template<class ...> class Predicate,
	 class ... Args2>
struct find_if_quinary_pred;

template<typename T1, typename T2, typename T3, typename T4,
	 template<class ...> class Predicate>
struct find_if_quinary_pred<T1, T2, T3, T4, Predicate>
  : std::integral_constant<std::size_t, 0>
{};

template<typename T1, typename T2, typename T3, typename T4,
	 template<class ...T> class Predicate,
	 class Head, class ... Tail>
struct find_if_quinary_pred<T1, T2, T3, T4, Predicate, Head, Tail...>
  : std::conditional <
  Predicate<Head, T1,T2,T3,T4>::type::value,
  std::integral_constant<std::size_t, 0>,
  std::integral_constant <
    std::size_t, 1 + find_if_quinary_pred<T1,T2,T3,T4, Predicate, Tail...>::type::value
    >
  >::type
{};

template <typename T1, typename T2, typename T3, typename T4,
	  template <class... T> class Predicate,
	  class... Args>
using find_if_quinary_pred_t = typename find_if_quinary_pred<T1, T2, T3, T4,
							     Predicate,
							     Args...>::type;

}}} // namespace

#endif  // MPL_VARIADIC_FIND_IF_QUINARY_PRED_HPP_
