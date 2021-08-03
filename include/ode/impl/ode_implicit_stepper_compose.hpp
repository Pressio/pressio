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

template<typename tag, typename ...Args>
struct ImplicitCompose;

////////////////////////////////////////
/// CRANK-NICOLSON stepper
////////////////////////////////////////
// 1. ImplicitCompose<tag, StateType, res_t, jac_t, system_t>;
// 2. ImplicitCompose<tag, StateType, res_t, jac_t, system_t, void, res_pol, jac_pol>;
template<
  typename StateType,
  typename ResidualType,
  typename JacobianType,
  typename SystemType
  >
struct ImplicitCompose<
  ::pressio::ode::implicitmethods::CrankNicolson,
  mpl::enable_if_t<
    ::pressio::ode::implicit_state<StateType>::value and
    ::pressio::ode::implicit_residual<ResidualType>::value and
    ::pressio::ode::implicit_jacobian<JacobianType>::value and
    ::pressio::ode::continuous_time_system_with_user_provided_jacobian<SystemType>::value
  >,
  StateType, ResidualType, JacobianType, SystemType>
{
  using ScalarType = typename ::pressio::Traits<StateType>::scalar_type;

  static_assert(::pressio::are_scalar_compatible<
   StateType, ResidualType, JacobianType>::value,
   "state, residual and jacobian are not scalar compatible ");

  using ResidualPolicyType =
    ::pressio::ode::impl::ResidualStandardPolicyCrankNicolson<
    StateType, ResidualType>;
  using JacobianPolicyType =
    ::pressio::ode::impl::JacobianStandardPolicyCrankNicolson<
    StateType, JacobianType>;

  using type = StepperCrankNicolson<
    ScalarType, StateType, ResidualType, JacobianType,
    SystemType, ResidualPolicyType, JacobianPolicyType,
    true // policies are standard
    >;
};

// 2. when we pass a void before policies in place of aux stepper
template<
  typename StateType,
  typename ResidualType,
  typename JacobianType,
  typename SystemType,
  typename ResidualPolicyType,
  typename JacobianPolicyType>
struct ImplicitCompose<
  ::pressio::ode::implicitmethods::CrankNicolson,
  mpl::enable_if_t<
    ::pressio::ode::implicit_state<StateType>::value and
    ::pressio::ode::implicit_residual<ResidualType>::value and
    ::pressio::ode::implicit_jacobian<JacobianType>::value
    and
    ::pressio::ode::implicit_cranknicolson_residual_policy<
      ResidualPolicyType, StateType, ResidualType, SystemType,
      typename ::pressio::Traits<StateType>::scalar_type
      >::value
    and
    ::pressio::ode::implicit_cranknicolson_jacobian_policy<
      JacobianPolicyType, StateType, JacobianType, SystemType,
      typename ::pressio::Traits<StateType>::scalar_type
      >::value
  >,
  StateType, ResidualType, JacobianType,
  SystemType, void, ResidualPolicyType, JacobianPolicyType>
{
  using ScalarType = typename ::pressio::Traits<StateType>::scalar_type;

  static_assert(::pressio::are_scalar_compatible
   <StateType, ResidualType, JacobianType>::value,
   "state, residual and jacobian are not scalar compatible ");

  using type = StepperCrankNicolson<
    ScalarType, StateType, ResidualType, JacobianType,
    SystemType, const ResidualPolicyType &, const JacobianPolicyType &,
    false // policies are custom
    >;
};


////////////////////////////////////////
/// BDF1 stepper
////////////////////////////////////////

// potential cases
// 1. ImplicitCompose<bdf1, StateType, res_t, jac_t, system_t>;
// 2. ImplicitCompose<bdf1, StateType, res_t, jac_t, system_t, res_pol, jac_pol>;
// 3. ImplicitCompose<bdf1, StateType, res_t, jac_t, system_t, void, res_pol, jac_pol>;

// 1.
// when we have standard policies, the system must be checked for API
// because the standard policies call some specific methods on the system,
// so the system must meet some requirements
template<
  typename StateType,
  typename ResidualType,
  typename JacobianType,
  typename SystemType
  >
struct ImplicitCompose<
  ::pressio::ode::implicitmethods::BDF1,
  mpl::enable_if_t<
    ::pressio::ode::implicit_state<StateType>::value and
    ::pressio::ode::implicit_residual<ResidualType>::value and
    ::pressio::ode::implicit_jacobian<JacobianType>::value and
    ::pressio::ode::continuous_time_system_with_user_provided_jacobian<SystemType>::value
  >,
  StateType, ResidualType, JacobianType, SystemType>
{
  using ScalarType = typename ::pressio::Traits<StateType>::scalar_type;
  static_assert(::pressio::are_scalar_compatible<
   StateType, ResidualType, JacobianType>::value,
   "state, residual and jacobian are not scalar compatible ");

  using ResidualPolicyType =
    ::pressio::ode::impl::ResidualStandardPolicyBdf<
    StateType, ResidualType>;
  using JacobianPolicyType =
    ::pressio::ode::impl::JacobianStandardPolicyBdf<
    StateType, JacobianType>;

  using type = StepperBDF1<
    ScalarType, StateType, ResidualType, JacobianType,
    SystemType, ResidualPolicyType, JacobianPolicyType,
    true // policies are standard
    >;
};

// 2.
// when we have custom policies, we do not check the system, because in general,
// one could not use it since the stepper only knows about the policies,
// and these can themselves contain the system object.
// in this case, we dont need to checck if system has certain requirements.
// it is user's responsibility to make sure the policies
// are doing the right thing.
template<
  typename StateType,
  typename ResidualType,
  typename JacobianType,
  typename SystemType,
  typename ResidualPolicyType,
  typename JacobianPolicyType
  >
struct ImplicitCompose<
  ::pressio::ode::implicitmethods::BDF1,
  mpl::enable_if_t<
    ::pressio::ode::implicit_state<StateType>::value and
    ::pressio::ode::implicit_residual<ResidualType>::value and
    ::pressio::ode::implicit_jacobian<JacobianType>::value
    and
    ::pressio::ode::implicit_euler_residual_policy<
      ResidualPolicyType, StateType, ResidualType, SystemType,
      typename ::pressio::Traits<StateType>::scalar_type>::value
    and
    ::pressio::ode::implicit_euler_jacobian_policy<
      JacobianPolicyType, StateType, JacobianType, SystemType,
      typename ::pressio::Traits<StateType>::scalar_type>::value
  >,
  StateType, ResidualType, JacobianType,
  SystemType, ResidualPolicyType, JacobianPolicyType>
{
  using ScalarType = typename ::pressio::Traits<StateType>::scalar_type;

  static_assert(::pressio::are_scalar_compatible
   <StateType, ResidualType, JacobianType>::value,
   "state, residual and jacobian are not scalar compatible ");

  using type = StepperBDF1<
    ScalarType, StateType, ResidualType, JacobianType,
    SystemType, const ResidualPolicyType &, const JacobianPolicyType&,
    false // policies are custom
    >;
};

// 3. when we pass a void before policies in place of aux stepper
template<
  typename StateType,
  typename ResidualType,
  typename JacobianType,
  typename SystemType,
  typename ResidualPolicyType,
  typename JacobianPolicyType>
struct ImplicitCompose<
  ::pressio::ode::implicitmethods::BDF1,
  mpl::enable_if_t<
    ::pressio::ode::implicit_state<StateType>::value and
    ::pressio::ode::implicit_residual<ResidualType>::value and
    ::pressio::ode::implicit_jacobian<JacobianType>::value
    and
    ::pressio::ode::implicit_euler_residual_policy<
      ResidualPolicyType, StateType, ResidualType, SystemType,
      typename ::pressio::Traits<StateType>::scalar_type
      >::value
    and
    ::pressio::ode::implicit_euler_jacobian_policy<
      JacobianPolicyType, StateType, JacobianType, SystemType,
      typename ::pressio::Traits<StateType>::scalar_type
      >::value
  >,
  StateType, ResidualType, JacobianType,
  SystemType, void, ResidualPolicyType, JacobianPolicyType>
{
  using ScalarType = typename ::pressio::Traits<StateType>::scalar_type;

  static_assert(::pressio::are_scalar_compatible
   <StateType, ResidualType, JacobianType>::value,
   "state, residual and jacobian are not scalar compatible ");

  using type = StepperBDF1<
    ScalarType, StateType, ResidualType, JacobianType,
    SystemType, const ResidualPolicyType &, const JacobianPolicyType &,
    false // policies are custom
    >;
};


////////////////////////////////////////
/// BDF2 stepper
////////////////////////////////////////

// ImplicitCompose<bdf2, StateType, res_t, jac_t, system_t, aux_stepper_t>;
// ImplicitCompose<bdf2, StateType, res_t, jac_t, system_t, aux_stepper_t, res_pol, jac_pol>;
template<
  typename StateType,
  typename ResidualType,
  typename JacobianType,
  typename SystemType,
  typename aux_stepper_t
  >
struct ImplicitCompose<
  ::pressio::ode::implicitmethods::BDF2,
  mpl::enable_if_t<
    ::pressio::ode::implicit_state<StateType>::value and
    ::pressio::ode::implicit_residual<ResidualType>::value and
    ::pressio::ode::implicit_jacobian<JacobianType>::value and
    ::pressio::ode::continuous_time_system_with_user_provided_jacobian<SystemType>::value and
    ::pressio::ode::auxiliary_stepper_for_bdf2<aux_stepper_t>::value
  >,
  StateType, ResidualType, JacobianType, SystemType, aux_stepper_t>
{
  using ScalarType = typename ::pressio::Traits<StateType>::scalar_type;
  static_assert(::pressio::are_scalar_compatible<
   StateType, ResidualType, JacobianType>::value,
   "state, residual and jacobian are not scalar compatible ");

  using ResidualPolicyType =
    ::pressio::ode::impl::ResidualStandardPolicyBdf<
    StateType, ResidualType>;
  using JacobianPolicyType =
    ::pressio::ode::impl::JacobianStandardPolicyBdf<
    StateType, JacobianType>;

  using type = StepperBDF2<
    ScalarType, StateType, ResidualType, JacobianType,
    SystemType, aux_stepper_t,
    ResidualPolicyType, JacobianPolicyType,
    true // policies are standard
    >;
};

template<
  typename StateType,
  typename ResidualType,
  typename JacobianType,
  typename SystemType,
  typename aux_stepper_t,
  typename ResidualPolicyType,
  typename JacobianPolicyType
  >
struct ImplicitCompose<
  ::pressio::ode::implicitmethods::BDF2,
  mpl::enable_if_t<
    ::pressio::ode::implicit_state<StateType>::value and
    ::pressio::ode::implicit_residual<ResidualType>::value and
    ::pressio::ode::implicit_jacobian<JacobianType>::value and
    ::pressio::ode::auxiliary_stepper_for_bdf2<aux_stepper_t>::value
    and
    ::pressio::ode::implicit_bdf2_residual_policy<
      ResidualPolicyType, StateType, ResidualType, SystemType,
      typename ::pressio::Traits<StateType>::scalar_type>::value
    and
    ::pressio::ode::implicit_bdf2_jacobian_policy<
      JacobianPolicyType, StateType, JacobianType, SystemType,
      typename ::pressio::Traits<StateType>::scalar_type>::value
  >,
  StateType, ResidualType, JacobianType,
  SystemType, aux_stepper_t,
  ResidualPolicyType, JacobianPolicyType>
{
  using ScalarType = typename ::pressio::Traits<StateType>::scalar_type;
  static_assert(::pressio::are_scalar_compatible
   <StateType, ResidualType, JacobianType>::value,
   "state, residual and jacobian are not scalar compatible ");

  using type = StepperBDF2<
    ScalarType, StateType, ResidualType, JacobianType,
    SystemType, aux_stepper_t,
    const ResidualPolicyType &, const JacobianPolicyType &,
    false // policies are custom
    >;
};


////////////////////////////////////////
/// ImplicitCompose Arbitrary stepper
////////////////////////////////////////

// ImplicitCompose< StateType, res_t, jac_t, system_t, OrderSetter, TotNumStatesSetter>;
// if we are here, it means we only need auxiliary states and not auxiliary rhs
template<
  typename StateType,
  typename ResidualType,
  typename JacobianType,
  typename SystemType,
  typename order_setter_t,
  typename tot_n_setter_t
  >
struct ImplicitCompose<
  ::pressio::ode::implicitmethods::Arbitrary,
  mpl::enable_if_t<
    ::pressio::ode::implicit_state<StateType>::value and
    ::pressio::ode::implicit_residual<ResidualType>::value and
    ::pressio::ode::implicit_jacobian<JacobianType>::value and
    ::pressio::ode::discrete_time_system_with_user_provided_jacobian<SystemType>::value and
    ::pressio::ode::IsStepperOrderSetter<order_setter_t>::value and
    ::pressio::ode::IsStepperTotalNumStatesSetter<tot_n_setter_t>::value
  >,
  StateType, ResidualType, JacobianType, SystemType, order_setter_t, tot_n_setter_t>
{
  using ScalarType = typename ::pressio::Traits<StateType>::scalar_type;
  static_assert(::pressio::are_scalar_compatible<
   StateType, ResidualType, JacobianType>::value,
   "state, residual and jacobian are not scalar compatible ");

  using ResidualPolicyType =
    ::pressio::ode::impl::ResidualStandardDiscreteTimePolicy<
    StateType, ResidualType>;
  using JacobianPolicyType =
    ::pressio::ode::impl::JacobianStandardDiscreteTimePolicy<
    StateType, JacobianType>;

  using type = StepperArbitrary<
    ScalarType, StateType, ResidualType, JacobianType,
    SystemType, order_setter_t, tot_n_setter_t,
    ResidualPolicyType, JacobianPolicyType,
    true // policies are standard
    >;
};

// ImplicitCompose< StateType, res_t, jac_t, system_t, OrderSetter, TotNumStatesSetter, res_pol, jac_pol>;
// if we are here, it means we only need auxiliary states and not auxiliary rhs
template<
  typename StateType,
  typename ResidualType,
  typename JacobianType,
  typename SystemType,
  typename order_setter_t,
  typename tot_n_setter_t,
  typename ResidualPolicyType,
  typename JacobianPolicyType
>
struct ImplicitCompose<
  ::pressio::ode::implicitmethods::Arbitrary,
  mpl::enable_if_t<
    ::pressio::ode::implicit_state<StateType>::value and
    ::pressio::ode::implicit_residual<ResidualType>::value and
    ::pressio::ode::implicit_jacobian<JacobianType>::value and
    ::pressio::ode::IsStepperOrderSetter<order_setter_t>::value and
    ::pressio::ode::IsStepperTotalNumStatesSetter<tot_n_setter_t>::value
    and
    ::pressio::ode::implicit_residual_policy<
      ResidualPolicyType, ::pressio::ode::implicitmethods::Arbitrary,
      tot_n_setter_t::value - 1, // number of wanted auxiliary states
      0, // number of wanted auxiliary rhs (this is because of the template specialization)
      StateType, ResidualType, SystemType,
      typename ::pressio::Traits<StateType>::scalar_type>::value
    and
    ::pressio::ode::implicit_jacobian_policy<
      JacobianPolicyType, ::pressio::ode::implicitmethods::Arbitrary,
      tot_n_setter_t::value - 1, // number of wanted auxiliary states
      StateType, JacobianType, SystemType,
      typename ::pressio::Traits<StateType>::scalar_type>::value
  >,
  StateType, ResidualType, JacobianType, SystemType, order_setter_t, tot_n_setter_t,
  ResidualPolicyType, JacobianPolicyType>
{
  using ScalarType = typename ::pressio::Traits<StateType>::scalar_type;
  static_assert(::pressio::are_scalar_compatible<
   StateType, ResidualType, JacobianType>::value,
   "state, residual and jacobian are not scalar compatible ");

  using type = StepperArbitrary<
    ScalarType, StateType, ResidualType, JacobianType,
    SystemType, order_setter_t, tot_n_setter_t,
    const ResidualPolicyType &, const JacobianPolicyType &,
    false // policies are custom
    >;
};

}}}
#endif  // ODE_IMPLICIT_IMPL_ODE_IMPLICIT_STEPPER_ImplicitCompose_IMPL_HPP_
