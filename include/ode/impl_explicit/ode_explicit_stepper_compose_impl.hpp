/*
//@HEADER
// ************************************************************************
//
// ode_explicit_stepper_compose_impl.hpp
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

#ifndef ODE_EXPLICIT_IMPL_ODE_EXPLICIT_STEPPER_COMPOSE_IMPL_HPP_
#define ODE_EXPLICIT_IMPL_ODE_EXPLICIT_STEPPER_COMPOSE_IMPL_HPP_

#include "ode_explicit_euler_stepper_impl.hpp"
#include "ode_explicit_runge_kutta4_stepper_impl.hpp"
#include "ode_explicit_adams_bashforth2_stepper_impl.hpp"

namespace pressio{ namespace ode{ namespace explicitmethods{ namespace impl{

template<typename tag>
struct ImplSelector
{
  template<bool, typename ...Args> using type = void;
};

template<>
struct ImplSelector<::pressio::ode::explicitmethods::Euler>
{
  template<bool is_standard_policy, typename ...Args>
  using type = ::pressio::ode::explicitmethods::impl::ExplicitEulerStepper
    <Args..., is_standard_policy>;
};

template<>
struct ImplSelector<::pressio::ode::explicitmethods::RungeKutta4>
{
  template<bool is_standard_policy, typename ...Args>
  using type = ::pressio::ode::explicitmethods::impl::ExplicitRungeKutta4Stepper
    <Args..., is_standard_policy>;
};

template<>
struct ImplSelector<::pressio::ode::explicitmethods::AdamsBashforth2>
{
  template<bool is_standard_policy, typename ...Args>
  using type = ::pressio::ode::explicitmethods::impl::ExplicitAdamsBashforth2Stepper
    <Args..., is_standard_policy>;
};


// --------------------------------------------------------------------------
// tag, state_type, system_t
// tag, state_type, system_t, policy
// --------------------------------------------------------------------------
template<
  typename tag,
  typename state_type,
  typename system_t,
  typename ...Args
  >
struct compose
{
  static_assert
  (::pressio::ode::continuous_time_system_with_at_least_velocity<system_t>::value,
   "The system passed to the ExplicitStepper does not meet the required API");

  static_assert
  (::pressio::ode::explicit_state<state_type>::value,
   "Invalid state type for explicit time stepping");

  using scalar_t = typename state_type::traits::scalar_t;
  using state_t	= state_type;
  using velocity_t = state_type;

  // typedef for standard velocity policy
  // (just typedef, only used if the user does not pass a user-defined policy)
  using standard_velocity_policy_t =
    ::pressio::ode::explicitmethods::policy::VelocityStandardPolicy<state_t>;

  // check Args if a user-defined velocity policy is passed
  using ic3 = ::pressio::mpl::variadic::find_if_quinary_pred_t<
    scalar_t, state_type, velocity_t, system_t,
    ::pressio::ode::explicit_velocity_policy, Args...>;
  using velocity_policy_t =
    ::pressio::mpl::variadic::at_or_t<standard_velocity_policy_t, ic3::value, Args...>;
  static constexpr bool is_standard_policy = std::is_same
    <standard_velocity_policy_t, velocity_policy_t>::value;

  // implementation class type
  using type =
    mpl::conditional_t<
    std::is_same<standard_velocity_policy_t, velocity_policy_t>::value,
    typename ImplSelector<tag>::template type<
      is_standard_policy, scalar_t, state_t, system_t,
      velocity_t, velocity_policy_t
      >,
    typename ImplSelector<tag>::template type<
      is_standard_policy, scalar_t, state_t, system_t,
      velocity_t, const velocity_policy_t &
      >
    >;
};

}}}} // end namespace pressio::ode::explicitmethods::impl
#endif  // ODE_EXPLICIT_IMPL_ODE_EXPLICIT_STEPPER_COMPOSE_IMPL_HPP_
