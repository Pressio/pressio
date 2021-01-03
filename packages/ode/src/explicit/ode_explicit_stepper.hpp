/*
//@HEADER
// ************************************************************************
//
// ode_explicit_stepper.hpp
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

#ifndef ODE_EXPLICIT_ODE_EXPLICIT_STEPPER_HPP_
#define ODE_EXPLICIT_ODE_EXPLICIT_STEPPER_HPP_

#include "./impl/ode_explicit_stepper_compose_impl.hpp"

namespace pressio{ namespace ode{

template<typename stepper_tag, typename ...Args>
using ExplicitStepper =
  typename ::pressio::ode::explicitmethods::impl::compose<
  stepper_tag,
  typename std::remove_cv<typename std::remove_reference<Args>::type>::type...
  >::type;

template<typename state_type, typename system_type, typename ...Args>
auto createForwardEulerStepper(const state_type & state,
			       const system_type & system,
			       Args && ...args)
  -> ExplicitStepper<explicitmethods::Euler, state_type,
		     system_type, Args...>
{
  using type = ExplicitStepper
    <explicitmethods::Euler, state_type, system_type, Args...>;
  return type(state, system, std::forward<Args>(args)...);
};

template<typename state_type, typename system_type, typename ...Args>
auto createRungeKutta4Stepper(const state_type & state,
			      const system_type & system,
			      Args && ...args)
  -> ExplicitStepper<explicitmethods::RungeKutta4, state_type,
		     system_type, Args...>
{
  using type = ExplicitStepper
    <explicitmethods::RungeKutta4, state_type, system_type, Args...>;
  return type(state, system, std::forward<Args>(args)...);
};

template<typename state_type, typename system_type, typename ...Args>
auto createAdamsBashforth2Stepper(const state_type & state,
			      const system_type & system,
			      Args && ...args)
  -> ExplicitStepper<explicitmethods::AdamsBashforth2, state_type,
		     system_type, Args...>
{
  using type = ExplicitStepper
    <explicitmethods::AdamsBashforth2, state_type, system_type, Args...>;
  return type(state, system, std::forward<Args>(args)...);
};

}} // end namespace pressio::ode
#endif  // ODE_EXPLICIT_ODE_EXPLICIT_STEPPER_HPP_
