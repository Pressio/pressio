/*
//@HEADER
// ************************************************************************
//
// ode_find_if_legitimate_implicit_jacobian_policy.hpp
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

#ifndef ODE_POLICIES_META_FIND_IF_LEGITIMATE_IMPLICIT_JACOBIAN_POLICY_HPP_
#define ODE_POLICIES_META_FIND_IF_LEGITIMATE_IMPLICIT_JACOBIAN_POLICY_HPP_

namespace pressio{ namespace ode{ namespace meta {

template<
  typename name_tag,
  typename state_t,
  typename jacobian_t,
  typename system_t,
  typename scalar_t,
  class ... Args2
  >
struct find_if_legitimate_implicit_jacobian_policy;

template<
  typename name_tag,
  typename state_t,
  typename jacobian_t,
  typename system_t,
  typename scalar_t
  >
struct find_if_legitimate_implicit_jacobian_policy<
  name_tag, state_t, jacobian_t, system_t, scalar_t
  > : std::integral_constant<std::size_t, 0>{};


template<
  typename name_tag,
  typename state_t,
  typename jacobian_t,
  typename system_t,
  typename scalar_t,
  class Head, class ... Tail
  >
struct find_if_legitimate_implicit_jacobian_policy<
  name_tag, state_t, jacobian_t, system_t, scalar_t,
  Head, Tail...
  >
  : std::conditional <
  is_legitimate_implicit_jacobian_policy
  <Head, name_tag,
   state_t, jacobian_t, system_t, scalar_t
   >::type::value,
  std::integral_constant<std::size_t, 0>,
  std::integral_constant <
    std::size_t, 1 +
    find_if_legitimate_implicit_jacobian_policy
    <name_tag, state_t, jacobian_t, system_t, scalar_t, Tail...>::type::value
    >
  >::type
{};

template <typename name_tag,
	  typename state_t,
	  typename jacobian_t,
	  typename system_t,
	  typename scalar_t,
	  class... Args>
using find_if_legitimate_implicit_jacobian_policy_t =
  typename find_if_legitimate_implicit_jacobian_policy<
  name_tag, state_t, jacobian_t, system_t, scalar_t, Args...>::type;


}}} // namespace pressio::ode::meta
#endif
