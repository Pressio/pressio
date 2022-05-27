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

template<class SystemType>
struct ExplicitCompose
{
  using sys_type = mpl::remove_cvref_t<SystemType>;

  // these static asserts are here because it allows to produce
  // readable error messages to the user if something is wrong.
  // If we were to us sfinae here, we would get a huge mess
  // for errors, so it is better to leave these are static asserts.
  // Ideally a system should be a concept. We cannot do that yet.

  static_assert
  (::pressio::ode::continuous_time_system_with_at_least_velocity<sys_type>::value,
   "To create an explicit stepper your system class MUST meet at least the velocity API");

  // if we get here, it means the above static assertion passed so
  // the system meets at least the velocity API
  using state_type = typename sys_type::state_type;
  static_assert(::pressio::ode::explicit_state<state_type>::value,
		"Invalid state type for explicit stepper");

  using velocity_type = typename sys_type::velocity_type;
  static_assert(::pressio::ode::explicit_velocity<velocity_type>::value,
		"Invalid velocity type for explicit time stepping");

  // it is very important to use "SystemType" as template arg
  // because that it the right type carrying how we store the system
  // since SystemType comes from the deduction below
  using mass_matrix_type = typename find_mass_matrix_if_any_or_noop<sys_type>::type;
  using type = ExplicitStepper<state_type, SystemType, velocity_type, mass_matrix_type>;
};

template<class SystemType>
auto create_explicit_stepper(StepScheme name, SystemType && system)
{

  // When user passes a non-temporary system object, SystemType is
  // deduced to be a reference, so the concrete stepper class
  // composed inside the ExplicitCompose will be composed such
  // that it will hold a reference to the provided system arg.
  // When the user passes system to be a temporary object,
  // SystemType will be deduced so that the stepper will
  // hold an **instance** of the system that
  // is move-constructed (if applicable) from the system argument.
  // So here it is important to use SystemType as template argument
  // for ExplicitCompose and NOT SystemType &&.
  using ReturnType = typename impl::ExplicitCompose<SystemType>::type;

  if (name == StepScheme::ForwardEuler){
    return ReturnType(ode::ForwardEuler(), std::forward<SystemType>(system));
  }

  else if (name == StepScheme::RungeKutta4){
    return ReturnType(ode::RungeKutta4(), std::forward<SystemType>(system));
  }

  else if (name == StepScheme::AdamsBashforth2){
    return ReturnType(ode::AdamsBashforth2(), std::forward<SystemType>(system));
  }

  else if (name == StepScheme::SSPRungeKutta3){
    return ReturnType(ode::SSPRungeKutta3(), std::forward<SystemType>(system));
  }

  else{
    throw std::runtime_error("ode:: create_explicit_stepper: invalid StepScheme enum value");
  }

}

}}}
#endif  // ODE_STEPPERS_IMPL_ODE_EXPLICIT_COMPOSE_HPP_
