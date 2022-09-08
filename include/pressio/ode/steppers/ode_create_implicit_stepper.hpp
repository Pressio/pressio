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
  Since we cannot yet use c++20 concepts, constraints should be enforced
  via e.g. SFINAE but that would yield bad error messages sometimes.
  So for now we decide to use, where possible, static asserts to have readable error messages.
*/

template<
  class SystemType,
  mpl::enable_if_t<
    ::pressio::ode::SystemWithRhsAndJacobian<SystemType>::value, int > = 0
  >
auto create_implicit_stepper(StepScheme name,
			     const SystemType & system)
{
  constexpr bool customPolicyCase = false;
  return impl::create_implicit_stepper_impl<customPolicyCase>(name, system);
}

template<
  class SystemType,
  class MassMatrixOperatorType,
  mpl::enable_if_t<
        ::pressio::ode::SystemWithRhsAndJacobian<SystemType>::value
    && (::pressio::ode::MassMatrixOperator<mpl::remove_cvref_t<MassMatrixOperatorType>>::value
     || ::pressio::ode::ConstantMassMatrixOperator<mpl::remove_cvref_t<MassMatrixOperatorType>>::value),
    int > = 0
  >
auto create_implicit_stepper(StepScheme name,
			     const SystemType & system,
			     MassMatrixOperatorType && massMatrixOperator)
{
  constexpr bool customPolicyCase = false;
  return impl::create_implicit_stepper_impl<
    customPolicyCase>(name, system, std::forward<MassMatrixOperatorType>(massMatrixOperator));
}

template<
  class ResidualJacobianPolicyType,
  mpl::enable_if_t<
    ::pressio::ode::ImplicitResidualJacobianPolicy<
      mpl::remove_cvref_t<ResidualJacobianPolicyType>>::value, int
    > = 0
  >
auto create_implicit_stepper(StepScheme name,
			     ResidualJacobianPolicyType && rjPolicy)
{
  constexpr bool customPolicyCase = true;
  return impl::create_implicit_stepper_impl<customPolicyCase>
    (name, std::forward<ResidualJacobianPolicyType>(rjPolicy));
}

// num of states as template arg constructs the arbitrary stepper
template<int TotalNumberOfDesiredStates, class SystemType>
auto create_implicit_stepper(SystemType && system)
{
  // the following should be a constraint
  using sys_type = mpl::remove_cvref_t<SystemType>;
  static_assert(::pressio::ode::FullyDiscreteSystemWithJacobian<
		sys_type, TotalNumberOfDesiredStates>::value,
		"The system passed does not meet the FullyDiscrete API");

  using sys_type = mpl::remove_cvref_t<SystemType>;
  using independent_variable_type = typename sys_type::independent_variable_type;
  using state_type = typename sys_type::state_type;
  using residual_type = typename sys_type::discrete_residual_type;
  using jacobian_type = typename sys_type::discrete_jacobian_type;

  using stepper_type = impl::StepperArbitrary<
    TotalNumberOfDesiredStates, SystemType, independent_variable_type,
    state_type, residual_type, jacobian_type
    >;

  return stepper_type(std::forward<SystemType>(system));
}

//
// auxiliary API
//
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
