/*
//@HEADER
// ************************************************************************
//
// ode_implicit_stepper_compose_impl.hpp
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

#ifndef ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_COMPOSE_IMPL_HPP_
#define ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_COMPOSE_IMPL_HPP_

#include "ode_implicit_discrete_time_residual.hpp"
#include "ode_implicit_discrete_time_jacobian.hpp"
#include "ode_implicit_policy_residual.hpp"
#include "ode_implicit_policy_jacobian.hpp"
#include "ode_implicit_stepper_euler.hpp"
#include "ode_implicit_stepper_bdf2.hpp"
#include "ode_implicit_stepper_cranknicolson.hpp"
#include "ode_implicit_stepper_arbitrary.hpp"

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

template<class ... Args>
struct ImplicitComposeAssertValidPolicies;

template<
  class ResidualPolicyType,
  class JacobianPolicyType,
  class ScalarType,
  class StateType
  >
struct ImplicitComposeAssertValidPolicies<
  ::pressio::ode::BDF1, ResidualPolicyType, JacobianPolicyType, ScalarType, StateType
  >
{
  static_assert(::pressio::ode::implicit_euler_residual_policy<
		mpl::remove_cvref_t<ResidualPolicyType>, StateType, ScalarType>::value,
		"Invalid residual policy for BDF1");
  static_assert(::pressio::ode::implicit_euler_jacobian_policy<
		mpl::remove_cvref_t<JacobianPolicyType>, StateType, ScalarType>::value,
		"Invalid jacobian policy for BDF1");
  static constexpr bool value = true;
};

template<
  class ResidualPolicyType,
  class JacobianPolicyType,
  class ScalarType,
  class StateType
  >
struct ImplicitComposeAssertValidPolicies<
  ::pressio::ode::BDF2, ResidualPolicyType, JacobianPolicyType, ScalarType, StateType
  >
{
  static_assert(::pressio::ode::implicit_bdf2_residual_policy<
		mpl::remove_cvref_t<ResidualPolicyType>, StateType, ScalarType>::value,
		"Invalid residual policy for BDF2");
  static_assert(::pressio::ode::implicit_bdf2_jacobian_policy<
		mpl::remove_cvref_t<JacobianPolicyType>, StateType, ScalarType>::value,
		"Invalid jacobian policy for BDF2");
  static constexpr bool value = true;
};

template<
  class ResidualPolicyType,
  class JacobianPolicyType,
  class ScalarType,
  class StateType
  >
struct ImplicitComposeAssertValidPolicies<
  ::pressio::ode::CrankNicolson, ResidualPolicyType, JacobianPolicyType, ScalarType, StateType
  >
{
  static_assert(::pressio::ode::implicit_cranknicolson_residual_policy<
		mpl::remove_cvref_t<ResidualPolicyType>, StateType, ScalarType>::value,
		"Invalid residual policy for CrankNicolson");
  static_assert(::pressio::ode::implicit_cranknicolson_jacobian_policy<
		mpl::remove_cvref_t<JacobianPolicyType>, StateType, ScalarType>::value,
		"Invalid jacobian policy for CrankNicolson");
  static constexpr bool value = true;
};


template<class Tag> struct ImplicitComposeConcreteStepper;

template<> struct ImplicitComposeConcreteStepper<::pressio::ode::BDF1>{
  template<bool b, class ...Args> using type = StepperBDF1<Args..., b>;
};

template<> struct ImplicitComposeConcreteStepper<::pressio::ode::BDF2>{
  template<bool b, class ...Args> using type = StepperBDF2<Args..., b>;
};

template<>
struct ImplicitComposeConcreteStepper<::pressio::ode::CrankNicolson>{
  template<bool b, class ...Args> using type = StepperCrankNicolson<Args..., b>;
};

template<class Tag> struct ImplicitComposeConcretePolicies;

template<>
struct ImplicitComposeConcretePolicies<::pressio::ode::BDF1>{
  template<class ...Args> using residual_policy_type = ResidualStandardPolicyBdf<Args...>;
  template<class ...Args> using jacobian_policy_type = JacobianStandardPolicy<Args...>;
};

template<>
struct ImplicitComposeConcretePolicies<::pressio::ode::BDF2>{
  template<class ...Args> using residual_policy_type = ResidualStandardPolicyBdf<Args...>;
  template<class ...Args> using jacobian_policy_type = JacobianStandardPolicy<Args...>;
};

template<>
struct ImplicitComposeConcretePolicies<::pressio::ode::CrankNicolson>{
  template<class ...Args> using residual_policy_type = ResidualStandardPolicyCrankNicolson<Args...>;
  template<class ...Args> using jacobian_policy_type = JacobianStandardPolicy<Args...>;
};


template<class ...Args>
struct ImplicitCompose{
  using type = void;
};

////////////////////////////////////////
/// BDF1 or BDF2 stepper
////////////////////////////////////////

// 0. ImplicitCompose<tag, enable, SystemType, StateType>
template<class TagType, class SystemType, class StateType>
struct ImplicitCompose<
  TagType,
  mpl::enable_if_t<
       std::is_same<TagType, ::pressio::ode::BDF1>::value
    or std::is_same<TagType, ::pressio::ode::BDF2>::value
    or std::is_same<TagType, ::pressio::ode::CrankNicolson>::value
    >,
  SystemType, StateType>
{
  static_assert( std::is_lvalue_reference<SystemType>::value,
  "The system instance must be an lvalue reference");

  static_assert
  (::pressio::ode::continuous_time_system_with_user_provided_jacobian<mpl::remove_cvref_t<SystemType>>::value,
   "The system passed does not meet the required API");

  using ResidualType = typename mpl::remove_cvref_t<SystemType>::velocity_type;
  using JacobianType = typename mpl::remove_cvref_t<SystemType>::jacobian_type;
  static_assert
  (ImplicitComposeAssertValidStateResJac<StateType, ResidualType, JacobianType>::value, "");

  using ScalarType = typename ::pressio::Traits<StateType>::scalar_type;
  using ResidualPolicyType = typename ImplicitComposeConcretePolicies<TagType>::template residual_policy_type<SystemType, StateType, ResidualType>;
  using JacobianPolicyType = typename ImplicitComposeConcretePolicies<TagType>::template jacobian_policy_type<SystemType, StateType, JacobianType>;

  using type = typename ImplicitComposeConcreteStepper<TagType>::template  type<
    true, ScalarType, StateType, ResidualType, JacobianType, ResidualPolicyType, JacobianPolicyType>;
};

// 1. ImplicitCompose<tag, enable, StateType, ResPol,  JacPol>
template<
  class TagType,
  class StateType,
  class ResidualPolicyType,
  class JacobianPolicyType
  >
struct ImplicitCompose<
  TagType,
  mpl::enable_if_t<
       std::is_same<TagType, ::pressio::ode::BDF1>::value
    or std::is_same<TagType, ::pressio::ode::BDF2>::value
    or std::is_same<TagType, ::pressio::ode::CrankNicolson>::value
    >,
  StateType, ResidualPolicyType, JacobianPolicyType>
{
  using ScalarType = typename ::pressio::Traits<StateType>::scalar_type;
  static_assert(ImplicitComposeAssertValidPolicies<
    TagType, ResidualPolicyType, JacobianPolicyType, ScalarType, StateType>::value, "");

  using ResidualType = typename mpl::remove_cvref_t<ResidualPolicyType>::residual_type;
  using JacobianType = typename mpl::remove_cvref_t<JacobianPolicyType>::jacobian_type;
  static_assert
  (ImplicitComposeAssertValidStateResJac<StateType, ResidualType, JacobianType>::value, "");

  using type = typename ImplicitComposeConcreteStepper<TagType>::template  type<
    false, ScalarType, StateType, ResidualType, JacobianType, ResidualPolicyType, JacobianPolicyType>;
};


template<typename ...Args>
using ImplicitCompose_t = typename ImplicitCompose<Args...>::type;

template<
  class TagType,
  class SystemType,
  class StateType,
  class ReturnType = ImplicitCompose_t<TagType, void, SystemType, StateType>
  >
ReturnType create_stepper_impl(SystemType && system, const StateType & state){
  return ReturnType(state, std::forward<SystemType>(system));
};

template<
  class TagType,
  class StateType,
  class ResidualPolicyType,
  class JacobianPolicyType,
  class ReturnType = ImplicitCompose_t<TagType, void, StateType, ResidualPolicyType, JacobianPolicyType>
  >
ReturnType create_stepper_impl(const StateType & state, ResidualPolicyType && rPol, JacobianPolicyType && jPol)
{
  return ReturnType(state, std::forward<ResidualPolicyType>(rPol), std::forward<JacobianPolicyType>(jPol));
};

////////////////////////////////////////
/// ImplicitCompose Arbitrary stepper
////////////////////////////////////////
template<int num_states, class SystemType, class StateType>
struct ImplicitComposeArb
{
  static_assert
  (::pressio::ode::discrete_time_system_with_user_provided_jacobian<
   mpl::remove_cvref_t<SystemType>, num_states>::value,
   "The system passed does not meet the required API");

  static_assert
  (::pressio::ode::implicit_state<StateType>::value,
   "Invalid state type for implicit stepper");

  static_assert
  (std::is_same<StateType, typename mpl::remove_cvref_t<SystemType>::state_type>::value,
   "Incompatible StateType and state_type alias deduced from the system class");

  using ResidualType = typename mpl::remove_cvref_t<SystemType>::discrete_time_residual_type;
  using JacobianType = typename mpl::remove_cvref_t<SystemType>::discrete_time_jacobian_type;
  static_assert(::pressio::ode::implicit_residual<ResidualType>::value,
    "Invalid residual type for implicit time stepping");
  static_assert(::pressio::ode::implicit_jacobian<JacobianType>::value,
    "Invalid jacobian type for implicit time stepping");
  static_assert(::pressio::are_scalar_compatible<StateType, ResidualType, JacobianType>::value,
   "state, residual and jacobian are not scalar compatible ");

  using ScalarType = typename ::pressio::Traits<StateType>::scalar_type;
  using type = StepperArbitrary<
    num_states, ScalarType, StateType, ResidualType, JacobianType, SystemType
    >;
};


}}}
#endif  // ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_ImplicitCompose_IMPL_HPP_
