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

#include "ode_implicit_discrete_time_residual.hpp"
#include "ode_implicit_discrete_time_jacobian.hpp"
#include "ode_implicit_policy_residual.hpp"
#include "ode_implicit_policy_jacobian.hpp"
#include "ode_implicit_stepper_arbitrary.hpp"
#include "ode_implicit_stepper_standard.hpp"

namespace pressio{ namespace ode{ namespace impl{

template<class StateType, class ResidualType,class JacobianType>
struct ImplicitComposeAssertValidStateResJac
{
  static_assert(::pressio::ode::implicit_state<StateType>::value,
    "Invalid state type for implicit time stepping");
  static_assert(::pressio::ode::implicit_residual<ResidualType>::value,
    "Invalid residual type for implicit time stepping");
  static_assert(::pressio::ode::implicit_jacobian<JacobianType>::value,
    "Invalid jacobian type for implicit time stepping");
  static_assert(::pressio::are_scalar_compatible<StateType, ResidualType, JacobianType>::value,
   "state, residual and jacobian are not scalar compatible ");

  static constexpr bool value = true;
};

template<class ResidualPolicyType, class JacobianPolicyType>
struct ImplicitComposeAssertValidPolicies
{
  static_assert(::pressio::ode::implicit_residual_policy<
		ResidualPolicyType >::value, "Invalid residual policy");

  static_assert(::pressio::ode::implicit_jacobian_policy<
		JacobianPolicyType >::value, "Invalid jacobian policy");

  static_assert(std::is_same<
		typename ResidualPolicyType::time_type,
		typename JacobianPolicyType::time_type
		>::value,
		"Residual and jacobian policies have mismatching time_type");

  static_assert( std::is_same<
		 typename ResidualPolicyType::state_type,
		 typename JacobianPolicyType::state_type
		 >::value,
		 "Mismatching nested state_type alias in residual and jacobian policies");

  static constexpr bool value = true;
};

////////////////////////////////////////
/// Arbitrary stepper
////////////////////////////////////////
template<int num_states, class SystemType>
struct ImplicitComposeArb
{

  using sys_type = mpl::remove_cvref_t<SystemType>;
  static_assert
  (::pressio::ode::discrete_time_system_with_user_provided_jacobian<
   sys_type, num_states>::value,
   "The system passed does not meet the required API");

  using state_type = typename sys_type::state_type;
  static_assert(::pressio::ode::implicit_state<state_type>::value,
		"Invalid state type for implicit stepper");

  using residual_type = typename sys_type::discrete_time_residual_type;
  using jacobian_type = typename sys_type::discrete_time_jacobian_type;
  static_assert(::pressio::ode::implicit_residual<residual_type>::value,
		"Invalid residual type for implicit time stepping");
  static_assert(::pressio::ode::implicit_jacobian<jacobian_type>::value,
		"Invalid jacobian type for implicit time stepping");

  using time_type = typename sys_type::time_type;
  using type = StepperArbitrary<
    num_states, time_type, state_type, residual_type, jacobian_type, SystemType
    >;
};

////////////////////////////////////////
/// standard stepper
////////////////////////////////////////
template<class ...Args>
struct ImplicitCompose{
  using type = void;
};

template<class SystemType>
struct ImplicitCompose<SystemType>
{
  using sys_type = mpl::remove_cvref_t<SystemType>;
  static_assert(::pressio::ode::continuous_time_system_with_user_provided_jacobian<
		sys_type>::value, "The system passed does not meet the required API");

  using time_type     = typename sys_type::time_type;
  using state_type    = typename sys_type::state_type;
  using residual_type = typename sys_type::velocity_type;
  using jacobian_type = typename sys_type::jacobian_type;
  static_assert(ImplicitComposeAssertValidStateResJac<
		state_type, residual_type, jacobian_type>::value, "");

  using mass_matrix_type = typename find_mass_matrix_if_any_or_noop<sys_type>::type;

  using residual_policy_type = ResidualStandardPolicy<
    SystemType, time_type, state_type, residual_type, mass_matrix_type>;
  using jacobian_policy_type = JacobianStandardPolicy<
    SystemType, time_type, state_type, jacobian_type, mass_matrix_type>;

  using type = StepperRt<time_type, state_type, residual_type, jacobian_type,
			 residual_policy_type, jacobian_policy_type, true>;
};

template<class ResidualPolicyType, class JacobianPolicyType>
struct ImplicitCompose<ResidualPolicyType, JacobianPolicyType>
{
  using residual_policy_type = typename mpl::remove_cvref_t<ResidualPolicyType>;
  using jacobian_policy_type = typename mpl::remove_cvref_t<JacobianPolicyType>;
  static_assert(ImplicitComposeAssertValidPolicies<
		residual_policy_type, jacobian_policy_type>::value, "");
  // time and state types are same for residual and jacobian policies
  using time_type  = typename residual_policy_type::time_type;
  using state_type = typename residual_policy_type::state_type;

  using residual_type = typename residual_policy_type::residual_type;
  using jacobian_type = typename jacobian_policy_type::jacobian_type;
  static_assert(ImplicitComposeAssertValidStateResJac<
		state_type, residual_type, jacobian_type>::value, "");

  using type = StepperRt<time_type, state_type, residual_type, jacobian_type,
			 ResidualPolicyType, jacobian_policy_type, false>;
};

template<class SystemType>
auto create_implicit_stepper_impl(StepScheme name, SystemType && system)
{
  using return_t = typename ImplicitCompose<SystemType>::type;

  if (name == StepScheme::BDF1){
    return return_t(::pressio::ode::BDF1(), system);
  }

  else if (name == StepScheme::BDF2){
    return return_t(::pressio::ode::BDF2(), system);
  }

  else if (name == StepScheme::CrankNicolson){
    return return_t(::pressio::ode::CrankNicolson(), system);
  }

  else{
    throw std::runtime_error("ode:: create_implicit_stepper: invalid StepScheme enum value");
  }
}

template<class ResidualPolicyType, class JacobianPolicyType>
auto create_implicit_stepper_impl(StepScheme name,
				  ResidualPolicyType && rPol,
				  JacobianPolicyType && jPol)
{
  using return_t = typename ImplicitCompose<ResidualPolicyType,
					    JacobianPolicyType>::type;

  if (name == StepScheme::BDF1){
    return return_t(::pressio::ode::BDF1(),
		      std::forward<ResidualPolicyType>(rPol),
		      std::forward<JacobianPolicyType>(jPol));
  }

  else if (name == StepScheme::BDF2){
    return return_t(::pressio::ode::BDF2(),
		      std::forward<ResidualPolicyType>(rPol),
		      std::forward<JacobianPolicyType>(jPol));
  }

  else if (name == StepScheme::CrankNicolson){
    return return_t(::pressio::ode::CrankNicolson(),
		      std::forward<ResidualPolicyType>(rPol),
		      std::forward<JacobianPolicyType>(jPol));
  }

  else{
    throw std::runtime_error("ode:: create_implicit_stepper: invalid StepScheme enum value");
  }
}

}}}
#endif  // ODE_STEPPERS_IMPL_ODE_IMPLICIT_STEPPER_COMPOSE_HPP_
