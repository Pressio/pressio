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
#include "ode_explicit_stepper_with_mass_matrix.hpp"
#include "ode_trivial_mass_matrix.hpp"

namespace pressio{ namespace ode{ namespace impl{

template<class SystemType, class = void>
struct ExplicitStepperImplClass{
  template<class ...Args> using type = ExplicitStepperNoMassMatrixImpl<Args...>;
};

template<class SystemType>
struct ExplicitStepperImplClass<
  SystemType, mpl::enable_if_t<system_has_constant_mass_matrix_api<SystemType>::value>
  >
{
  template<class ...Args> using type = ExplicitStepperWithMassMatrixImpl<
    true /*signal that we have mass matrix but it is constant */,
    typename SystemType::mass_matrix_type, Args...>;
};

template<class SystemType>
struct ExplicitStepperImplClass<
  SystemType, mpl::enable_if_t<system_has_potentially_varying_mass_matrix_api<SystemType>::value>
  >
{
  template<class ...Args> using type = ExplicitStepperWithMassMatrixImpl<
    false /*signal that we have mass matrix but it is NOT const */,
    typename SystemType::mass_matrix_type, Args...>;
};

template<class SystemType>
struct ExplicitCompose
{
  using sys_type = mpl::remove_cvref_t<SystemType>;

  // these static asserts are here because it allows to produce
  // readable error messages to the user if something is wrong.
  // If we were to us sfinae here, we would get a huge mess
  // for errors, so it is better to leave these are static asserts.
  // Ideally we should use concepts for better errors. We cannot do that yet.

  static_assert
  (   ::pressio::ode::OdeSystemWithRhs<sys_type>::value
   || ::pressio::ode::OdeSystemWithRhsAndMassMatrix<sys_type>::value
   || ::pressio::ode::OdeSystemWithRhsAndJacobian<sys_type>::value
   || ::pressio::ode::OdeSystemWithRhsAndConstantMassMatrix<sys_type>::value
   || ::pressio::ode::OdeSystemComplete<sys_type>::value
   || ::pressio::ode::OdeSystemCompleteWithConstantMassMatrix<sys_type>::value,
   "explicit stepper: your system class does not meet any valid concept");

  // if we get here, it means the above static assertion passed so
  // the system meets at least the right_hand_side API
  using state_type = typename sys_type::state_type;
  using right_hand_side_type = typename sys_type::right_hand_side_type;

  // it is very important to use "SystemType" as template arg
  // because that it the right type carrying how we store the system
  // since SystemType comes from the deduction below
  using ind_var_type = typename sys_type::independent_variable_type;
  using type = typename ExplicitStepperImplClass<sys_type>::template type<
    state_type, ind_var_type, SystemType, right_hand_side_type>;
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
