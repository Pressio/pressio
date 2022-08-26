/*
//@HEADER
// ************************************************************************
//
// ode_implicit_stepper_compose.hpp
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

#ifndef ODE_STEPPERS_IMPL_ODE_IMPLICIT_STEPPER_COMPOSE_HPP_
#define ODE_STEPPERS_IMPL_ODE_IMPLICIT_STEPPER_COMPOSE_HPP_

#include "ode_implicit_discrete_residual.hpp"
#include "ode_implicit_discrete_jacobian.hpp"

#include "ode_implicit_policy_residual_jacobian.hpp"
#include "ode_implicit_policy_residual_jacobian_with_mass_matrix.hpp"

#include "ode_implicit_stepper_arbitrary.hpp"
#include "ode_implicit_stepper_standard.hpp"

namespace pressio{ namespace ode{ namespace impl{

////////////////////////////////////////
/// standard stepper
// ////////////////////////////////////////

template<bool usingCustomPoliciy, class ...Args>
struct ImplicitCompose{
  using type = void;
};

// partial specialize for NO custom policy, NO mass matrix
template<class SystemType>
struct ImplicitCompose<false, SystemType>
{
  using sys_type = std::decay_t<SystemType>;

  // if we get here, the system has AT LEAST rhs, jacobian
  using independent_variable_type = typename sys_type::independent_variable_type;
  using state_type    = typename sys_type::state_type;
  using residual_type = typename sys_type::right_hand_side_type;
  using jacobian_type = typename sys_type::jacobian_type;

  // we HAVE to pass SystemType here as template because it carries
  // the correct qualifiers or the policy instantiation won't work
  using rj_policy_type = ResidualJacobianStandardPolicy<
    SystemType, independent_variable_type, state_type, residual_type, jacobian_type>;
  // using jacobian_policy_type = JacobianStandardPolicy<
  //   SystemType, independent_variable_type, state_type, jacobian_type>;

  using type = ImplicitStepperStandardImpl<
    independent_variable_type, state_type, residual_type, jacobian_type,
    rj_policy_type,
    false, // false because there is no mass matrix
    true /*this is using default policies*/>;
};

// partial specialize for NO custom policy, WITH mass matrix
template<class SystemType, class MassMatrixOpType>
struct ImplicitCompose<false, SystemType, MassMatrixOpType>
{
  using sys_type = std::decay_t<SystemType>;

  // if we get here, the system has AT LEAST rhs, jacobian
  using independent_variable_type = typename sys_type::independent_variable_type;
  using state_type    = typename sys_type::state_type;
  using residual_type = typename sys_type::right_hand_side_type;
  using jacobian_type = typename sys_type::jacobian_type;

  using rj_policy_type = ResidualJacobianWithMassMatrixStandardPolicy<
    ConstantMassMatrixOperator<std::decay_t<MassMatrixOpType>>::value,
    SystemType, independent_variable_type,
    state_type, residual_type, jacobian_type,
    MassMatrixOpType>;

  using type = ImplicitStepperStandardImpl<
    independent_variable_type, state_type, residual_type, jacobian_type,
    rj_policy_type,
    true, /*with mass matrix*/
    true /*this is using default policies*/>;
};

// partial specialize for custom policy, NO mass matrix
template<class ResidualJacobianPolicyType>
struct ImplicitCompose<true, ResidualJacobianPolicyType>
{
  using rj_policy_type = std::decay_t<ResidualJacobianPolicyType>;

  using independent_variable_type  = typename rj_policy_type::independent_variable_type;
  using state_type    = typename rj_policy_type::state_type;
  using residual_type = typename rj_policy_type::residual_type;
  using jacobian_type = typename rj_policy_type::jacobian_type;

  using type = ImplicitStepperStandardImpl<
    independent_variable_type, state_type, residual_type, jacobian_type,
    ResidualJacobianPolicyType,
    false, /*no mass matrix*/
    false /*this is using custom policies*/>;
};

template<bool customPolicy, class ... Args>
auto create_implicit_stepper_impl(StepScheme name,
				  Args && ... args)
{

  using return_t = typename ImplicitCompose<customPolicy, Args...>::type;
  if (name == StepScheme::BDF1){
    return return_t(::pressio::ode::BDF1(), std::forward<Args>(args)...);
  }
  else if (name == StepScheme::BDF2){
    return return_t(::pressio::ode::BDF2(), std::forward<Args>(args)...);
  }
  else if (name == StepScheme::CrankNicolson){
    return return_t(::pressio::ode::CrankNicolson(), std::forward<Args>(args)...);
  }
  else{
    throw std::runtime_error("ode:: create_implicit_stepper: invalid StepScheme enum value");
  }
}

}}}
#endif  // ODE_STEPPERS_IMPL_ODE_IMPLICIT_STEPPER_COMPOSE_HPP_
