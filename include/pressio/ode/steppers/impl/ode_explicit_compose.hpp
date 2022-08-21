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

namespace pressio{ namespace ode{ namespace impl{

template<class SystemType, class ... Args>
struct ExplicitComposeReturnType;

template<class SystemType>
struct ExplicitComposeReturnType<void, SystemType>
{
  using sys_type = std::decay_t<SystemType>;
  using ind_var_type = typename sys_type::independent_variable_type;
  using state_type = typename sys_type::state_type;
  using right_hand_side_type = typename sys_type::right_hand_side_type;

  // it is very important to use "SystemType" as template arg
  // because that it the right type carrying how we store the system
  // since SystemType comes from the deduction below
  using type = ExplicitStepperNoMassMatrixImpl<
    state_type, ind_var_type, SystemType, right_hand_side_type>;
};

template<class SystemType, class MassMatrixOpType>
struct ExplicitComposeReturnType<
  void, SystemType, MassMatrixOpType
  >
{
  using sys_type = std::decay_t<SystemType>;
  using state_type = typename sys_type::state_type;
  using right_hand_side_type = typename sys_type::right_hand_side_type;
  using ind_var_type = typename sys_type::independent_variable_type;

  // it is very important to use "SystemType" MassMatrixOpType as template arg
  // because that it the right type carrying how we store the system
  // since SystemType comes from the deduction below
  using type = ExplicitStepperWithMassMatrixImpl<
    ConstantMassMatrixOperator<std::decay_t<MassMatrixOpType>>::value,
    MassMatrixOpType, state_type, ind_var_type, SystemType, right_hand_side_type>;
};

template<class SystemType, class ... Args>
auto create_explicit_stepper(StepScheme name,
			     SystemType && system,
			     Args && ... args)
{

  // Use SystemType as template argument for ExplicitComposeReturnType
  // and NOT SystemType && for the following reason:
  // when user passes a non-temporary system object, SystemType is
  // deduced to be a reference, so the concrete stepper class
  // composed inside the ExplicitComposeReturnType will be composed such
  // that it will hold a reference to the provided system arg.
  // When the user passes system to be a temporary object,
  // SystemType will be deduced so that the stepper will
  // hold an **instance** of the system that
  // is move-constructed (if applicable) from the system argument.
  using ReturnType = typename impl::ExplicitComposeReturnType<
    void, SystemType, Args...>::type;

  if (name == StepScheme::ForwardEuler){
    return ReturnType(ode::ForwardEuler(),
		      std::forward<SystemType>(system),
		      std::forward<Args>(args)...);
  }

  else if (name == StepScheme::RungeKutta4){
    return ReturnType(ode::RungeKutta4(),
		      std::forward<SystemType>(system),
		      std::forward<Args>(args)...);
  }

  else if (name == StepScheme::AdamsBashforth2){
    return ReturnType(ode::AdamsBashforth2(),
		      std::forward<SystemType>(system),
		      std::forward<Args>(args)...);
  }

  else if (name == StepScheme::SSPRungeKutta3){
    return ReturnType(ode::SSPRungeKutta3(),
		      std::forward<SystemType>(system),
		      std::forward<Args>(args)...);
  }

  else{
    throw std::runtime_error("ode:: create_explicit_stepper: invalid StepScheme enum value");
  }

}

}}}
#endif  // ODE_STEPPERS_IMPL_ODE_EXPLICIT_COMPOSE_HPP_
