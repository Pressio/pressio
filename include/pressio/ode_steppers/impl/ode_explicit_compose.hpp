/*
//@HEADER
// ************************************************************************
//
// ode_explicit_compose.hpp
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

#ifndef ODE_STEPPERS_IMPL_ODE_EXPLICIT_COMPOSE_HPP_
#define ODE_STEPPERS_IMPL_ODE_EXPLICIT_COMPOSE_HPP_

#include "ode_explicit_stepper.hpp"

namespace pressio{ namespace ode{ namespace impl{

template<class StateType, class SystemType>
struct ExplicitCompose
{
  static_assert
  (::pressio::ode::continuous_time_system_with_at_least_velocity<mpl::remove_cvref_t<SystemType>>::value,
   "The system passed to the ExplicitStepper does not meet the required API");

  static_assert(::pressio::ode::explicit_state<StateType>::value,
		"Invalid state type for explicit stepper");
  static_assert
  (std::is_same<StateType, typename mpl::remove_cvref_t<SystemType>::state_type>::value,
   "Incompatible StateType and state_type alias deduced from the system class");

  using scalar_type   = typename ::pressio::Traits<StateType>::scalar_type;
  using scalar_type_from_system = typename mpl::remove_cvref_t<SystemType>::scalar_type;
  static_assert(std::is_same<scalar_type, scalar_type_from_system>::value,
		"State and system have inconsistent scalar types");

  using velocity_type = typename mpl::remove_cvref_t<SystemType>::velocity_type;
  static_assert(::pressio::ode::explicit_velocity<velocity_type>::value,
		"Invalid velocity type for explicit time stepping");

  // it is very important that the stepper is given "SystemType" because
  // that contains the right logic for how we store the system
  using type = ExplicitStepper<scalar_type, StateType, SystemType, velocity_type>;
};

template<class StateType, class SystemType>
auto create_explicit_stepper(StepScheme name,
			     const StateType & state,
			     SystemType && system)
{

  // note that here it is important to use SystemType and NOT SystemType &&.
  // When user passes a non-temporary system, SystemType is deduced to be a reference,
  // so the concrete stepper class composed inside the ExplicitCompose will
  // be composed such that it will hold a reference to the provided system arg.
  // When the user passes a temporary system object, SystemType will be correctly
  // deduced so that the stepper will hold an instance of the system that
  // is move-constructed (if applicable) from the system argument.
  using ReturnType = typename impl::ExplicitCompose<StateType, SystemType>::type;

  if (name == StepScheme::ForwardEuler){
    return ReturnType(ode::ForwardEuler(), state, std::forward<SystemType>(system));
  }

  else if (name == StepScheme::RungeKutta4){
    return ReturnType(ode::RungeKutta4(), state, std::forward<SystemType>(system));
  }

  else if (name == StepScheme::AdamsBashforth2){
    return ReturnType(ode::AdamsBashforth2(), state, std::forward<SystemType>(system));
  }

  else if (name == StepScheme::SSPRungeKutta3){
    return ReturnType(ode::SSPRungeKutta3(), state, std::forward<SystemType>(system));
  }

  else{
    throw std::runtime_error("ode:: create_explicit_stepper: invalid StepScheme enum value");
  }

};

}}}
#endif  // ODE_STEPPERS_IMPL_ODE_EXPLICIT_COMPOSE_HPP_
