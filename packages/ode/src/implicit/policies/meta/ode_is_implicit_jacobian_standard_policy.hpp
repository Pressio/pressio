/*
//@HEADER
// ************************************************************************
//
// ode_is_implicit_jacobian_standard_policy.hpp
//                     		      Pressio 
// Copyright 2019 National Technology & Engineering Solutions of Sandia,LLC 
//							      (NTESS)
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

#ifndef ODE_POLICIES_META_IS_IMPLICIT_JACOBIAN_STD_POLICY_HPP_
#define ODE_POLICIES_META_IS_IMPLICIT_JACOBIAN_STD_POLICY_HPP_

#include "../base/ode_jacobian_policy_base.hpp"
#include "../standard/ode_implicit_jacobian_standard_policy.hpp"

namespace pressio{ namespace ode{ namespace meta {

template<
  ::pressio::ode::ImplicitEnum whichone,
  typename policy_t,
  typename enable = void
  >
struct is_implicit_jacobian_standard_policy : std::false_type{};

template <
  template <typename...> class policy_t,
  typename... Args
  >
struct is_implicit_jacobian_standard_policy<
  ::pressio::ode::ImplicitEnum::Euler,
  policy_t<Args...>,
  typename std::enable_if<
    std::is_same<
      policy_t<Args...>,
      ode::policy::ImplicitJacobianStandardPolicy<
	Args...
	>
      >::value
    >::type > : std::true_type{};


template <
  template <typename...> class policy_t,
  typename... Args
  >
struct is_implicit_jacobian_standard_policy<
  ::pressio::ode::ImplicitEnum::BDF2,
  policy_t<Args...>,
  typename std::enable_if<
    std::is_same<
      policy_t<Args...>,
      ode::policy::ImplicitJacobianStandardPolicy<
	Args...
	>
      >::value
    >::type > : std::true_type{};


template<template <typename...> class policy_t, typename... Args>
using is_implicit_euler_jacobian_standard_policy =
  typename is_implicit_jacobian_standard_policy<
  ::pressio::ode::ImplicitEnum::Euler,policy_t<Args...>>::type;

template<template <typename...> class policy_t, typename... Args>
using is_implicit_bdf2_jacobian_standard_policy =
  typename is_implicit_jacobian_standard_policy<
  ::pressio::ode::ImplicitEnum::BDF2,policy_t<Args...>>::type;


}}} // namespace pressio::ode::meta
#endif
