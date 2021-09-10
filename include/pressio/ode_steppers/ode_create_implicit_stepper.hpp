/*
//@HEADER
// ************************************************************************
//
// ode_implicit_stepper.hpp
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

#ifndef ODE_IMPLICIT_ODE_IMPLICIT_STEPPER_HPP_
#define ODE_IMPLICIT_ODE_IMPLICIT_STEPPER_HPP_

#include "./impl/ode_implicit_stepper_compose.hpp"

namespace pressio{ namespace ode{

template<class StateType, class SystemType>
auto create_implicit_stepper(StepScheme name,
			     const StateType & state,
			     const SystemType & system)
  ->  decltype(impl::create_implicit_stepper_impl(name, system, state))
{
  return impl::create_implicit_stepper_impl(name, system, state);
};

template<class StateType, class ResidualPolicyType, class JacobianPolicyType>
auto create_implicit_stepper(StepScheme name,
			     const StateType & state,
			     ResidualPolicyType && rPol,
			     JacobianPolicyType && jPol)
  ->  decltype(impl::create_implicit_stepper_impl(name, state,
						  std::forward<ResidualPolicyType>(rPol),
						  std::forward<JacobianPolicyType>(jPol)))
{
  return impl::create_implicit_stepper_impl(name, state,
					    std::forward<ResidualPolicyType>(rPol),
					    std::forward<JacobianPolicyType>(jPol));
};

template<class ...Args>
auto create_bdf1_stepper(Args && ... args)
  -> decltype(create_implicit_stepper(StepScheme::BDF1, std::forward<Args>(args)...))
{
  return create_implicit_stepper(StepScheme::BDF1, std::forward<Args>(args)...);
};

template<class ...Args>
auto create_bdf2_stepper(Args && ... args)
  -> decltype(create_implicit_stepper(StepScheme::BDF2, std::forward<Args>(args)...))
{
  return create_implicit_stepper(StepScheme::BDF2, std::forward<Args>(args)...);
};

template<class ...Args>
auto create_cranknicolson_stepper(Args && ... args)
  -> decltype(create_implicit_stepper(StepScheme::CrankNicolson, std::forward<Args>(args)...))
{
  return create_implicit_stepper(StepScheme::CrankNicolson, std::forward<Args>(args)...);
};

// Arbitrary
template<
  int num_states,
  class StateType,
  class SystemType,
  class ReturnType = typename impl::ImplicitComposeArb<num_states, SystemType, StateType>::type
  >
ReturnType create_arbitrary_stepper(const StateType & state, SystemType && system)
{
  return ReturnType(state, std::forward<SystemType>(system));
};

}} // end namespace pressio::ode
#endif  // ODE_IMPLICIT_ODE_IMPLICIT_STEPPER_HPP_
