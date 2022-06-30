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

#include "ode_implicit_policy_residual.hpp"
#include "ode_implicit_policy_jacobian.hpp"
#include "ode_implicit_policy_residual_mass_matrix.hpp"
#include "ode_implicit_policy_jacobian_mass_matrix.hpp"

#include "ode_trivial_mass_matrix.hpp"
#include "ode_implicit_stepper_arbitrary.hpp"
#include "ode_implicit_stepper_standard.hpp"

namespace pressio{ namespace ode{ namespace impl{

template<class ResidualPolicyType, class JacobianPolicyType>
struct ImplicitComposeAssertValidPolicies
{
  static_assert(::pressio::ode::ImplicitResidualPolicy<
		ResidualPolicyType >::value, "Invalid residual policy");

  static_assert(::pressio::ode::ImplicitJacobianPolicy<
		JacobianPolicyType >::value, "Invalid jacobian policy");

  static_assert(std::is_same<
		typename ResidualPolicyType::independent_variable_type,
		typename JacobianPolicyType::independent_variable_type
		>::value,
		"Residual and jacobian policies have mismatching independent_variable_type");

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
  static_assert(::pressio::ode::FullyDiscreteSystemWithJacobian<
		sys_type, num_states>::value,
		"The system passed does not meet the required API");

  using state_type = typename sys_type::state_type;
  using residual_type = typename sys_type::discrete_residual_type;
  using jacobian_type = typename sys_type::discrete_jacobian_type;

  using independent_variable_type = typename sys_type::independent_variable_type;
  using type = StepperArbitrary<
    num_states, independent_variable_type, state_type, residual_type, jacobian_type, SystemType
    >;
};

////////////////////////////////////////
/// standard stepper
////////////////////////////////////////

template<class T, class = void>
struct PoliciesTypes{
  template<class ...Args>
  using r_policy_type = ResidualStandardPolicy<T, Args ...>;

  template<class ...Args>
  using j_policy_type = JacobianStandardPolicy<T, Args...>;
};

template<class T>
struct PoliciesTypes<
  T, mpl::enable_if_t<system_has_constant_mass_matrix_api< mpl::remove_cvref_t<T> >::value>
  >{
  using mass_mat_type = typename mpl::remove_cvref_t<T>::mass_matrix_type;
  template<class ...Args> using r_policy_type = ResidualWithMassMatrixStandardPolicy<
    true /*const mass matrix */, T, Args ..., mass_mat_type>;
  template<class ...Args> using j_policy_type = JacobianWithMassMatrixStandardPolicy<
    T, Args..., mass_mat_type>;
};

template<class T>
struct PoliciesTypes<
  T, mpl::enable_if_t<system_has_potentially_varying_mass_matrix_api< mpl::remove_cvref_t<T> >::value>
  >{
  using mass_mat_type = typename mpl::remove_cvref_t<T>::mass_matrix_type;
  template<class ...Args> using r_policy_type = ResidualWithMassMatrixStandardPolicy<
    false /*mass matrix is not const */, T, Args ..., mass_mat_type>;
  template<class ...Args> using j_policy_type = JacobianWithMassMatrixStandardPolicy<
    T, Args..., mass_mat_type>;
};

template<class ...Args>
struct ImplicitCompose{
  using type = void;
};

template<class SystemType>
struct ImplicitCompose<SystemType>
{
  using sys_type = mpl::remove_cvref_t<SystemType>;

  static_assert
  (   ::pressio::ode::SemiDiscreteSystemWithRhsAndJacobian<sys_type>::value
   || ::pressio::ode::CompleteSemiDiscreteSystem<sys_type>::value
   || ::pressio::ode::CompleteSemiDiscreteSystemWithConstantMassMatrix<sys_type>::value,
   "implicit stepper: your system class does not meet any valid concept");

  // if we get here, the system has AT LEAST rhs, jacobian
  using independent_variable_type = typename sys_type::independent_variable_type;
  using state_type    = typename sys_type::state_type;
  using residual_type = typename sys_type::right_hand_side_type;
  using jacobian_type = typename sys_type::jacobian_type;

  // we HAVE to pass SystemType here as template because it carries
  // the correct qualifiers or the policy instantiation won't work
  using residual_policy_type = typename PoliciesTypes<SystemType>::template r_policy_type<
    independent_variable_type, state_type, residual_type>;
  using jacobian_policy_type = typename PoliciesTypes<SystemType>::template j_policy_type<
    independent_variable_type, state_type, jacobian_type>;

  using type = ImplicitStepperStandardImpl<
    independent_variable_type, state_type, residual_type, jacobian_type,
    residual_policy_type, jacobian_policy_type,
    true /*this is using default policies*/>;
};

template<class ResidualPolicyType, class JacobianPolicyType>
struct ImplicitCompose<ResidualPolicyType, JacobianPolicyType>
{
  using residual_policy_type = typename mpl::remove_cvref_t<ResidualPolicyType>;
  using jacobian_policy_type = typename mpl::remove_cvref_t<JacobianPolicyType>;
  static_assert(ImplicitComposeAssertValidPolicies<
		residual_policy_type, jacobian_policy_type>::value, "");
  // ind_var and state types are same for residual and jacobian policies
  using independent_variable_type  = typename residual_policy_type::independent_variable_type;
  using state_type    = typename residual_policy_type::state_type;
  using residual_type = typename residual_policy_type::residual_type;
  using jacobian_type = typename jacobian_policy_type::jacobian_type;

  using type = ImplicitStepperStandardImpl<
    independent_variable_type, state_type, residual_type, jacobian_type,
    ResidualPolicyType, jacobian_policy_type,
    false /*this is using custom policies*/>;
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
  using return_t = typename ImplicitCompose<
    ResidualPolicyType, JacobianPolicyType>::type;

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
