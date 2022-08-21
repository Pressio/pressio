/*
//@HEADER
// ************************************************************************
//
// ode_create_implicit_stepper.hpp
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

#ifndef ODE_STEPPERS_ODE_CREATE_IMPLICIT_STEPPER_HPP_
#define ODE_STEPPERS_ODE_CREATE_IMPLICIT_STEPPER_HPP_

#include "./impl/ode_implicit_stepper_compose.hpp"

namespace pressio{ namespace ode{

/*
  below we use static asserts even for constraints but this is
  not fully correct because constraints should have an impact on the
  overload resolution read this:
    https://timsong-cpp.github.io/cppwp/n4861/structure#footnote-154

  Since we cannot yet use c++20 concepts, we should enforce these
  constraints via e.g. SFINAE but that would yield bad error messages.
  So for now we decide to use static asserts to have readable error messages.
  Another point that kind of justifies this here, for now, is that
  we have a simple enough overload set.
*/

template<class SystemType>
auto create_implicit_stepper(StepScheme name, const
			     SystemType & system)
{
  // the following should be a constraint
  static_assert
  (::pressio::ode::OdeRhsAndJacobianEvaluator<SystemType>::value,
   "implicit stepper: your system class does not meet any valid concept");

  return impl::create_implicit_stepper_implA(name, system);
}

template<
  class SystemType,
  class MassMatrixOperatorType,
  mpl::enable_if_t<
    ::pressio::ode::OdeRhsAndJacobianEvaluator<SystemType>::value
    && (::pressio::ode::MassMatrixOperator<std::decay_t<MassMatrixOperatorType>>::value
     || ::pressio::ode::ConstantMassMatrixOperator<std::decay_t<MassMatrixOperatorType>>::value),
    int > = 0
  >
auto create_implicit_stepper(StepScheme name,
			     const SystemType & system,
			     MassMatrixOperatorType && mmOperator)
{

  return impl::create_implicit_stepper_implA(name, system,
					    std::forward<MassMatrixOperatorType>(mmOperator));
}

template<
  class ResidualPolicyType,
  class JacobianPolicyType,
  mpl::enable_if_t<
    ::pressio::ode::ImplicitResidualPolicy<std::decay_t<ResidualPolicyType>>::value
    && ::pressio::ode::ImplicitJacobianPolicy<std::decay_t<JacobianPolicyType>>::value, int
    > = 0
  >
auto create_implicit_stepper(StepScheme name,
			     ResidualPolicyType && resPolicy,
			     JacobianPolicyType && jacPolicy)
{
  using residual_policy_type = std::decay_t<ResidualPolicyType>;
  using jacobian_policy_type = std::decay_t<JacobianPolicyType>;

  // the following two checks on nested typedefs are mandates
  // so it is fine to use static_assert
  static_assert(std::is_same<
		typename residual_policy_type::independent_variable_type,
		typename jacobian_policy_type::independent_variable_type
		>::value,
		"Residual and jacobian policies cannot have mismatching nested independent_variable_type");
  static_assert( std::is_same<
		 typename residual_policy_type::state_type,
		 typename jacobian_policy_type::state_type
		 >::value,
		"Residual and jacobian policies cannot have mismatching nested state_type");

  return impl::create_implicit_stepper_implB(name,
					    std::forward<ResidualPolicyType>(resPolicy),
					    std::forward<JacobianPolicyType>(jacPolicy));
}

// num of states as template arg constructs the arbitrary stepper
template<int num_states, class SystemType>
auto create_implicit_stepper(SystemType && system)
{
  // the following should be a constraint
  using sys_type = mpl::remove_cvref_t<SystemType>;
  static_assert(::pressio::ode::FullyDiscreteSystemWithJacobian<
		sys_type, num_states>::value,
		"The system passed does not meet the FullyDiscrete API");

  return typename impl::ImplicitComposeArb<
    num_states, SystemType>::type(std::forward<SystemType>(system));
}

template<class ...Args>
auto create_bdf1_stepper(Args && ... args){
  return create_implicit_stepper(StepScheme::BDF1,
				 std::forward<Args>(args)...);
}

template<class ...Args>
auto create_bdf2_stepper(Args && ... args){
  return create_implicit_stepper(StepScheme::BDF2,
				 std::forward<Args>(args)...);
}

template<class ...Args>
auto create_cranknicolson_stepper(Args && ... args){
  return create_implicit_stepper(StepScheme::CrankNicolson,
				 std::forward<Args>(args)...);
}

}} // end namespace pressio::ode
#endif  // ODE_STEPPERS_ODE_CREATE_IMPLICIT_STEPPER_HPP_
