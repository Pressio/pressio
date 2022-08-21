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

#include "ode_implicit_stepper_arbitrary.hpp"
#include "ode_implicit_stepper_standard.hpp"

namespace pressio{ namespace ode{ namespace impl{

////////////////////////////////////////
/// Arbitrary stepper
////////////////////////////////////////
template<int num_states, class SystemType>
struct ImplicitComposeArb
{
  using sys_type = mpl::remove_cvref_t<SystemType>;
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
// ////////////////////////////////////////

// template<class T, class = void>
// struct PoliciesTypes{
//   template<class ...Args>
//   using r_policy_type = ResidualStandardPolicy<T, Args ...>;

//   template<class ...Args>
//   using j_policy_type = JacobianStandardPolicy<T, Args...>;
// };

// template<class T>
// struct PoliciesTypes<
//   T, mpl::enable_if_t<system_has_constant_mass_matrix_api< mpl::remove_cvref_t<T> >::value>
//   >{
//   using mass_mat_type = typename mpl::remove_cvref_t<T>::mass_matrix_type;
//   template<class ...Args> using r_policy_type = ResidualWithMassMatrixStandardPolicy<
//     true /*const mass matrix */, T, Args ..., mass_mat_type>;
//   template<class ...Args> using j_policy_type = JacobianWithMassMatrixStandardPolicy<
//     T, Args..., mass_mat_type>;
// };

// template<class T>
// struct PoliciesTypes<
//   T, mpl::enable_if_t<system_has_potentially_varying_mass_matrix_api< mpl::remove_cvref_t<T> >::value>
//   >{
//   using mass_mat_type = typename mpl::remove_cvref_t<T>::mass_matrix_type;
//   template<class ...Args> using r_policy_type = ResidualWithMassMatrixStandardPolicy<
//     false /*mass matrix is not const */, T, Args ..., mass_mat_type>;
//   template<class ...Args> using j_policy_type = JacobianWithMassMatrixStandardPolicy<
//     T, Args..., mass_mat_type>;
// };

template<class ...Args>
struct ImplicitComposeA{
  using type = void;
};

template<class SystemType>
struct ImplicitComposeA<SystemType>
{
  using sys_type = std::decay_t<SystemType>;

  // if we get here, the system has AT LEAST rhs, jacobian
  using independent_variable_type = typename sys_type::independent_variable_type;
  using state_type    = typename sys_type::state_type;
  using residual_type = typename sys_type::right_hand_side_type;
  using jacobian_type = typename sys_type::jacobian_type;

  // we HAVE to pass SystemType here as template because it carries
  // the correct qualifiers or the policy instantiation won't work
  using residual_policy_type = ResidualStandardPolicy<
    SystemType, independent_variable_type, state_type, residual_type>;
  using jacobian_policy_type = JacobianStandardPolicy<
    SystemType, independent_variable_type, state_type, jacobian_type>;

  using type = ImplicitStepperStandardImpl<
    independent_variable_type, state_type, residual_type, jacobian_type,
    residual_policy_type, jacobian_policy_type,
    false,
    true /*this is using default policies*/>;
};

template<class SystemType, class MassMatrixOpType>
struct ImplicitComposeA<SystemType, MassMatrixOpType>
{
  using sys_type = std::decay_t<SystemType>;

  // if we get here, the system has AT LEAST rhs, jacobian
  using independent_variable_type = typename sys_type::independent_variable_type;
  using state_type    = typename sys_type::state_type;
  using residual_type = typename sys_type::right_hand_side_type;
  using jacobian_type = typename sys_type::jacobian_type;

  using residual_policy_type = ResidualWithMassMatrixStandardPolicy<
    ConstantMassMatrixOperator<std::decay_t<MassMatrixOpType>>::value,
    SystemType, independent_variable_type, state_type, residual_type, MassMatrixOpType>;

  using jacobian_policy_type = JacobianWithMassMatrixStandardPolicy<
    SystemType, independent_variable_type, state_type, jacobian_type,
    typename std::decay_t<MassMatrixOpType>::mass_matrix_type >;

  using type = ImplicitStepperStandardImpl<
    independent_variable_type, state_type, residual_type, jacobian_type,
    residual_policy_type, jacobian_policy_type,
    true, /*with mass matrix*/
    true /*this is using default policies*/>;
};

template<class ResidualPolicyType, class JacobianPolicyType>
struct ImplicitComposeB
{
  using residual_policy_type = typename mpl::remove_cvref_t<ResidualPolicyType>;
  using jacobian_policy_type = typename mpl::remove_cvref_t<JacobianPolicyType>;
  // ind_var and state types are same for residual and jacobian policies
  using independent_variable_type  = typename residual_policy_type::independent_variable_type;
  using state_type    = typename residual_policy_type::state_type;
  using residual_type = typename residual_policy_type::residual_type;
  using jacobian_type = typename jacobian_policy_type::jacobian_type;

  using type = ImplicitStepperStandardImpl<
    independent_variable_type, state_type, residual_type, jacobian_type,
    ResidualPolicyType, jacobian_policy_type,
    false, /*no mass matrix*/
    false /*this is using custom policies*/>;
};

template<class SystemType, class ... Args>
auto create_implicit_stepper_implA(StepScheme name,
				   SystemType && system,
				   Args && ... args)
{

  using return_t = typename ImplicitComposeA<SystemType, Args...>::type;
  if (name == StepScheme::BDF1){
    return return_t(::pressio::ode::BDF1(),
		    std::forward<SystemType>(system),
		    std::forward<Args>(args)...);
  }
  else if (name == StepScheme::BDF2){
    return return_t(::pressio::ode::BDF2(),
		    std::forward<SystemType>(system),
		    std::forward<Args>(args)...);
  }
  else if (name == StepScheme::CrankNicolson){
    return return_t(::pressio::ode::CrankNicolson(),
		    std::forward<SystemType>(system),
		    std::forward<Args>(args)...);
  }
  else{
    throw std::runtime_error("ode:: create_implicit_stepper: invalid StepScheme enum value");
  }
}

template<class ResidualPolicyType, class JacobianPolicyType>
auto create_implicit_stepper_implB(StepScheme name,
				   ResidualPolicyType && rPol,
				   JacobianPolicyType && jPol)
{
  using return_t = typename ImplicitComposeB<ResidualPolicyType, JacobianPolicyType>::type;

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
