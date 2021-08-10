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

#include "ode_implicit_policy_residual_bdf.hpp"
#include "ode_implicit_policy_jacobian_bdf.hpp"
#include "ode_implicit_policy_residual_crank_nicolson.hpp"
#include "ode_implicit_policy_jacobian_crank_nicolson.hpp"
#include "ode_implicit_policy_discrete_time_residual.hpp"
#include "ode_implicit_policy_discrete_time_jacobian.hpp"

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
  class StateType,
  class SystemType,
  class ResidualType,
  class JacobianType
  >
struct ImplicitComposeAssertValidPolicies<
  implicitmethods::BDF1,
  ResidualPolicyType, JacobianPolicyType,
  ScalarType, StateType, SystemType, ResidualType, JacobianType
  >
{
  static_assert(::pressio::ode::implicit_euler_residual_policy<
		mpl::remove_cvref_t<ResidualPolicyType>,
		StateType, ResidualType, SystemType, ScalarType>::value,
		"Invalid residual policy for BDF1");

  static_assert(::pressio::ode::implicit_euler_jacobian_policy<
		mpl::remove_cvref_t<JacobianPolicyType>,
		StateType, JacobianType, SystemType, ScalarType>::value,
		"Invalid jacobian policy for BDF1");
  static constexpr bool value = true;
};

template<
  class ResidualPolicyType,
  class JacobianPolicyType,
  class ScalarType,
  class StateType,
  class SystemType,
  class ResidualType,
  class JacobianType
  >
struct ImplicitComposeAssertValidPolicies<
  implicitmethods::BDF2,
  ResidualPolicyType, JacobianPolicyType,
  ScalarType, StateType, SystemType, ResidualType, JacobianType
  >
{
  static_assert(::pressio::ode::implicit_bdf2_residual_policy<
		mpl::remove_cvref_t<ResidualPolicyType>,
		StateType, ResidualType, SystemType, ScalarType>::value,
		"Invalid residual policy for BDF2");

  static_assert(::pressio::ode::implicit_bdf2_jacobian_policy<
		mpl::remove_cvref_t<JacobianPolicyType>,
		StateType, JacobianType, SystemType, ScalarType>::value,
		"Invalid jacobian policy for BDF2");
  static constexpr bool value = true;
};

template<
  class ResidualPolicyType,
  class JacobianPolicyType,
  class ScalarType,
  class StateType,
  class SystemType,
  class ResidualType,
  class JacobianType
  >
struct ImplicitComposeAssertValidPolicies<
  implicitmethods::CrankNicolson,
  ResidualPolicyType, JacobianPolicyType,
  ScalarType, StateType, SystemType, ResidualType, JacobianType
  >
{
  static_assert(::pressio::ode::implicit_cranknicolson_residual_policy<
		mpl::remove_cvref_t<ResidualPolicyType>,
		StateType, ResidualType, SystemType, ScalarType>::value,
		"Invalid residual policy for CrankNicolson");

  static_assert(::pressio::ode::implicit_cranknicolson_jacobian_policy<
		mpl::remove_cvref_t<JacobianPolicyType>,
		StateType, JacobianType, SystemType, ScalarType>::value,
		"Invalid jacobian policy for CrankNicolson");
  static constexpr bool value = true;
};


template<class Tag>
struct ImplicitComposeConcreteStepper;

template<>
struct ImplicitComposeConcreteStepper<implicitmethods::BDF1>
{
  template<bool b, class ...Args>
  using type = StepperBDF1<Args..., b>;
};

template<>
struct ImplicitComposeConcreteStepper<implicitmethods::BDF2>
{
  template<bool b, class ...Args>
  using type = StepperBDF2<Args..., b>;
};

template<>
struct ImplicitComposeConcreteStepper<implicitmethods::CrankNicolson>
{
  template<bool b, class ...Args>
  using type = StepperCrankNicolson<Args..., b>;
};


template<class Tag>
struct ImplicitComposeConcretePolicies;

template<>
struct ImplicitComposeConcretePolicies<implicitmethods::BDF1>{
  template<class ...Args> using residual_policy_type = ResidualStandardPolicyBdf<Args...>;
  template<class ...Args> using jacobian_policy_type = JacobianStandardPolicyBdf<Args...>;
};

template<>
struct ImplicitComposeConcretePolicies<implicitmethods::BDF2>{
  template<class ...Args> using residual_policy_type = ResidualStandardPolicyBdf<Args...>;
  template<class ...Args> using jacobian_policy_type = JacobianStandardPolicyBdf<Args...>;
};


template<>
struct ImplicitComposeConcretePolicies<implicitmethods::CrankNicolson>{
  template<class ...Args> using residual_policy_type = ResidualStandardPolicyCrankNicolson<Args...>;
  template<class ...Args> using jacobian_policy_type = JacobianStandardPolicyCrankNicolson<Args...>;
};


template<class ...Args>
struct ImplicitCompose{
  using type = void;
};

////////////////////////////////////////
/// BDF1 or BDF2 stepper
////////////////////////////////////////

// 0. ImplicitCompose<tag, enable, SystemType, StateType>
// 1. ImplicitCompose<tag, enable, SystemType, StateType, ResPol,  JacPol>
// 2. ImplicitCompose<tag, enable, SystemType, StateType, ResType, JacType, void>
// 3. ImplicitCompose<tag, enable, SystemType, StateType, ResType, JacType, ResPol, JacPol>

// 0.
template<class TagType, class SystemType, class StateType>
struct ImplicitCompose<
  TagType,
  mpl::enable_if_t<
    std::is_same<TagType, implicitmethods::BDF1>::value or
    std::is_same<TagType, implicitmethods::BDF2>::value or
    std::is_same<TagType, implicitmethods::CrankNicolson>::value
    >,
  SystemType, StateType>
{
  static_assert
  (::pressio::ode::continuous_time_system_with_user_provided_jacobian<SystemType>::value,
   "The system passed does not meet the required API");

  using ResidualType = typename SystemType::velocity_type;
  using JacobianType = typename SystemType::jacobian_type;
  static_assert
  (ImplicitComposeAssertValidStateResJac<StateType, ResidualType, JacobianType>::value, "");

  using ScalarType = typename ::pressio::Traits<StateType>::scalar_type;
  using ResidualPolicyType = typename ImplicitComposeConcretePolicies<TagType>::template residual_policy_type<StateType, ResidualType>;
  using JacobianPolicyType = typename ImplicitComposeConcretePolicies<TagType>::template jacobian_policy_type<StateType, JacobianType>;

  using type = typename ImplicitComposeConcreteStepper<TagType>::template  type<
    true, ScalarType, StateType, ResidualType, JacobianType,
    SystemType, ResidualPolicyType, JacobianPolicyType>;
};

// 1.
template<
  class TagType,
  class SystemType,
  class StateType,
  class ResidualPolicyType,
  class JacobianPolicyType
  >
struct ImplicitCompose<
  TagType,
  mpl::enable_if_t<
    std::is_same<TagType, implicitmethods::BDF1>::value or
    std::is_same<TagType, implicitmethods::BDF2>::value or
    std::is_same<TagType, implicitmethods::CrankNicolson>::value
    >,
  SystemType, StateType, ResidualPolicyType, JacobianPolicyType>
{

  static_assert
  (::pressio::ode::continuous_time_system_with_user_provided_jacobian<SystemType>::value,
   "The system passed does not meet the required API");

  using ResidualType = typename SystemType::velocity_type;
  using JacobianType = typename SystemType::jacobian_type;
  static_assert
  (ImplicitComposeAssertValidStateResJac<StateType, ResidualType, JacobianType>::value, "");
  using ScalarType = typename ::pressio::Traits<StateType>::scalar_type;

  static_assert(ImplicitComposeAssertValidPolicies<
		TagType, ResidualPolicyType, JacobianPolicyType,
		ScalarType, StateType, SystemType, ResidualType, JacobianType
		>::value, "");

  using type = typename ImplicitComposeConcreteStepper<TagType>::template  type<
    false, ScalarType, StateType, ResidualType, JacobianType,
    SystemType, ResidualPolicyType, JacobianPolicyType>;
};

// 2.
template<
  class TagType,
  class SystemType,
  class StateType,
  class ResidualType,
  class JacobianType
  >
struct ImplicitCompose<
  TagType,
  mpl::enable_if_t<
    std::is_same<TagType, implicitmethods::BDF1>::value or
    std::is_same<TagType, implicitmethods::BDF2>::value or
    std::is_same<TagType, implicitmethods::CrankNicolson>::value
    >,
  SystemType,
  StateType, ResidualType, JacobianType,
  void
  >
{

  static_assert
  (::pressio::ode::continuous_time_system_with_user_provided_jacobian<SystemType>::value,
   "The system passed does not meet the required API");
  static_assert
  (ImplicitComposeAssertValidStateResJac<StateType, ResidualType, JacobianType>::value, "");

  using ScalarType = typename ::pressio::Traits<StateType>::scalar_type;
  using ResidualPolicyType = typename ImplicitComposeConcretePolicies<TagType>::template residual_policy_type<StateType, ResidualType>;
  using JacobianPolicyType = typename ImplicitComposeConcretePolicies<TagType>::template jacobian_policy_type<StateType, JacobianType>;

  using type = typename ImplicitComposeConcreteStepper<TagType>::template type<
    true, ScalarType, StateType, ResidualType, JacobianType,
    SystemType, ResidualPolicyType, JacobianPolicyType>;
};

// 3.
// when we have custom policies and provided types for res and jac,
// we do not check the system, because in general,
// one could not use it since the stepper only knows about the policies,
// and these can themselves contain the system object.
// in this case, we dont need to check if system has certain requirements.
// it is user's responsibility to make sure the policies are doing the right thing.
template<
  class TagType,
  class SystemType,
  class StateType,
  class ResidualType,
  class JacobianType,
  class ResidualPolicyType,
  class JacobianPolicyType
  >
struct ImplicitCompose<
  TagType,
  mpl::enable_if_t<
    std::is_same<TagType, implicitmethods::BDF1>::value or
    std::is_same<TagType, implicitmethods::BDF2>::value or
    std::is_same<TagType, implicitmethods::CrankNicolson>::value
    >,
  SystemType,
  StateType, ResidualType, JacobianType,
  ResidualPolicyType, JacobianPolicyType
  >
{

  static_assert
  (ImplicitComposeAssertValidStateResJac<StateType, ResidualType, JacobianType>::value, "");
  using ScalarType = typename ::pressio::Traits<StateType>::scalar_type;

  static_assert(ImplicitComposeAssertValidPolicies<
		TagType, ResidualPolicyType, JacobianPolicyType,
		ScalarType, StateType, SystemType, ResidualType, JacobianType
		>::value, "");

  using type = typename ImplicitComposeConcreteStepper<TagType>::template type<
    false, ScalarType, StateType, ResidualType, JacobianType,
    SystemType, ResidualPolicyType, JacobianPolicyType>;
};


////////////////////////////////////////
/// ImplicitCompose Arbitrary stepper
////////////////////////////////////////

// 0. ImplicitCompose<tag, order, n_states, SystemType, StateType>
// 1. ImplicitCompose<tag, order, n_states, SystemType, StateType, ResPol,  JacPol>
// 2. ImplicitCompose<tag, order, n_states, SystemType, StateType, ResType, JacType, void>
// 3. ImplicitCompose<tag, order, n_states, SystemType, StateType, ResType, JacType, ResPol, JacPol>

// 0.
template<class Order, class NStates, class SystemType, class StateType>
struct ImplicitCompose<
  implicitmethods::Arbitrary, Order, NStates, SystemType, StateType>
{
  static_assert
  (::pressio::ode::discrete_time_system_with_user_provided_jacobian<SystemType>::value,
   "The system passed does not meet the required API");

  using ResidualType = typename SystemType::discrete_time_residual_type;
  using JacobianType = typename SystemType::discrete_time_jacobian_type;
  static_assert
  (ImplicitComposeAssertValidStateResJac<StateType, ResidualType, JacobianType>::value, "");

  using ScalarType = typename ::pressio::Traits<StateType>::scalar_type;
  using ResidualPolicyType = ResidualStandardDiscreteTimePolicy<StateType, ResidualType>;
  using JacobianPolicyType = JacobianStandardDiscreteTimePolicy<StateType, JacobianType>;

  using type = StepperArbitrary<
    Order::value, NStates::value,
    ScalarType, StateType, ResidualType, JacobianType,
    SystemType, ResidualPolicyType, JacobianPolicyType,
    true // policies are standard
    >;
};

// 1.
template<
  class Order,
  class NStates,
  class SystemType,
  class StateType,
  class ResidualPolicyType,
  class JacobianPolicyType
  >
struct ImplicitCompose<
  implicitmethods::Arbitrary, Order, NStates, SystemType, StateType, ResidualPolicyType, JacobianPolicyType
  >
{

  static_assert
  (::pressio::ode::discrete_time_system_with_user_provided_jacobian<SystemType>::value,
   "The system passed does not meet the required API");

  using ResidualType = typename SystemType::discrete_time_residual_type;
  using JacobianType = typename SystemType::discrete_time_jacobian_type;
  static_assert
  (ImplicitComposeAssertValidStateResJac<StateType, ResidualType, JacobianType>::value, "");
  using ScalarType = typename ::pressio::Traits<StateType>::scalar_type;

  static_assert(::pressio::ode::implicit_residual_policy<
		ResidualPolicyType,
		implicitmethods::Arbitrary,
		NStates::value - 1, // number of wanted auxiliary states
		0, // number of wanted auxiliary rhs (this is because of the template specialization)
		StateType, ResidualType, SystemType,
		ScalarType>::value,
		"Invalid residual policy for arbitrary stepper");

  static_assert(::pressio::ode::implicit_jacobian_policy<
		JacobianPolicyType,
		implicitmethods::Arbitrary,
		NStates::value - 1, // number of wanted auxiliary states
		StateType, JacobianType, SystemType,
		ScalarType>::value,
		"Invalid jacobian policy for arbitrary stepper");

  using type = StepperArbitrary<
    Order::value, NStates::value,
    ScalarType, StateType, ResidualType, JacobianType,
    SystemType, ResidualPolicyType, JacobianPolicyType,
    false
    >;
};

// 2.
template<
  class Order, class NStates,
  class SystemType,
  class StateType,
  class ResidualType,
  class JacobianType
  >
struct ImplicitCompose<
  implicitmethods::Arbitrary, Order, NStates, SystemType, StateType, ResidualType, JacobianType, void>
{
  static_assert
  (::pressio::ode::discrete_time_system_with_user_provided_jacobian<SystemType>::value,
   "The system passed does not meet the required API");
  static_assert
  (ImplicitComposeAssertValidStateResJac<StateType, ResidualType, JacobianType>::value, "");

  using ScalarType = typename ::pressio::Traits<StateType>::scalar_type;
  using ResidualPolicyType = ResidualStandardDiscreteTimePolicy<StateType, ResidualType>;
  using JacobianPolicyType = JacobianStandardDiscreteTimePolicy<StateType, JacobianType>;

  using type = StepperArbitrary<
    Order::value, NStates::value,
    ScalarType, StateType, ResidualType, JacobianType,
    SystemType, ResidualPolicyType, JacobianPolicyType,
    true // policies are standard
    >;
};


// 3.
template<
  class Order,
  class NStates,
  class SystemType,
  class StateType,
  class ResidualType,
  class JacobianType,
  class ResidualPolicyType,
  class JacobianPolicyType
  >
struct ImplicitCompose<
  implicitmethods::Arbitrary,
  Order, NStates,
  SystemType, StateType, ResidualType, JacobianType, ResidualPolicyType, JacobianPolicyType
  >
{
  static_assert
  (ImplicitComposeAssertValidStateResJac<StateType, ResidualType, JacobianType>::value, "");
  using ScalarType = typename ::pressio::Traits<StateType>::scalar_type;

  static_assert(::pressio::ode::implicit_residual_policy<
		ResidualPolicyType,
		implicitmethods::Arbitrary,
		NStates::value - 1, // number of wanted auxiliary states
		0, // number of wanted auxiliary rhs (this is because of the template specialization)
		StateType, ResidualType, SystemType,
		ScalarType>::value,
		"Invalid residual policy for arbitrary stepper");

  static_assert(::pressio::ode::implicit_jacobian_policy<
		JacobianPolicyType,
		implicitmethods::Arbitrary,
		NStates::value - 1, // number of wanted auxiliary states
		StateType, JacobianType, SystemType,
		ScalarType>::value,
		"Invalid jacobian policy for arbitrary stepper");

  using type = StepperArbitrary<
    Order::value, NStates::value,
    ScalarType, StateType, ResidualType, JacobianType,
    SystemType, ResidualPolicyType, JacobianPolicyType,
    false
    >;
};

template<class ...Args>
using ImplicitCompose_t = typename ImplicitCompose<Args...>::type;

template<
  class TagType,
  class SystemType,
  class StateType,
  class ReturnType = impl::ImplicitCompose_t<TagType, void, SystemType, StateType>
  >
ReturnType create_stepper_impl(const SystemType & system, const StateType & state){
  return ReturnType(state, system);
};

template<
  class TagType,
  class SystemType,
  class StateType,
  class ResidualPolicyType,
  class JacobianPolicyType,
  class ReturnType = impl::ImplicitCompose_t<TagType, void, SystemType, StateType, ResidualPolicyType, JacobianPolicyType>
  >
  ReturnType create_stepper_impl(const SystemType & system, const StateType & state, 
                                 ResidualPolicyType && rPol, JacobianPolicyType && jPol)
{
  return ReturnType(state, system,
        std::forward<ResidualPolicyType>(rPol),
        std::forward<JacobianPolicyType>(jPol));
};

template<
  class TagType,
  class ResidualType,
  class JacobianType,
  class SystemType,
  class StateType,
  class ReturnType = impl::ImplicitCompose_t<TagType, void, SystemType, StateType, ResidualType, JacobianType, void>
  >
ReturnType create_stepper_partial_deduction_impl(const SystemType & system, const StateType & state)
{
  return ReturnType(state, system);
};

template<
  class TagType,
  class ResidualType,
  class JacobianType,
  class SystemType,
  class StateType,
  class ResidualPolicyType,
  class JacobianPolicyType,
  class ReturnType = impl::ImplicitCompose_t<
    TagType, void, SystemType, StateType, ResidualType, JacobianType, ResidualPolicyType, JacobianPolicyType>
  >
ReturnType create_stepper_partial_deduction_impl(const SystemType & system, const StateType & state, 
                                                 ResidualPolicyType && rPol, JacobianPolicyType && jPol)
{
  return ReturnType(state, system,
        std::forward<ResidualPolicyType>(rPol),
        std::forward<JacobianPolicyType>(jPol));
};

}}}
#endif  // ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_ImplicitCompose_IMPL_HPP_
