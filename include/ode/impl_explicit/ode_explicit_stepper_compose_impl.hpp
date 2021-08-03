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
#include "ode_explicit_ssp_runge_kutta3_stepper_impl.hpp"

namespace pressio{ namespace ode{ namespace impl{

template<typename tag>
struct ExplicitImplSelector
{
  template<bool, typename ...Args> using type = void;
};

template<>
struct ExplicitImplSelector<::pressio::ode::explicitmethods::Euler>
{
  template<bool is_standard_policy, typename ...Args>
  using type = ::pressio::ode::impl::ExplicitEulerStepper
    <Args..., is_standard_policy>;
};

template<>
struct ExplicitImplSelector<::pressio::ode::explicitmethods::RungeKutta4>
{
  template<bool is_standard_policy, typename ...Args>
  using type = ::pressio::ode::impl::ExplicitRungeKutta4Stepper
    <Args..., is_standard_policy>;
};

template<>
struct ExplicitImplSelector<::pressio::ode::explicitmethods::AdamsBashforth2>
{
  template<bool is_standard_policy, typename ...Args>
  using type = ::pressio::ode::impl::ExplicitAdamsBashforth2Stepper
    <Args..., is_standard_policy>;
};

template<>
struct ExplicitImplSelector<::pressio::ode::explicitmethods::SSPRungeKutta3>
{
  template<bool is_standard_policy, typename ...Args>
  using type = ::pressio::ode::impl::ExplicitSSPRungeKutta3Stepper
    <Args..., is_standard_policy>;
};


template<class tag, class state_type, class system_t>
struct ExplicitComposeForDefaultPolicy
{
  static_assert
  (::pressio::ode::continuous_time_system_with_at_least_velocity<system_t>::value,
   "The system passed to the ExplicitStepper does not meet the required API");

  static_assert
  (::pressio::ode::explicit_state<state_type>::value,
   "Invalid state type for explicit time stepping");

  using scalar_type   = typename ::pressio::traits<state_type>::scalar_type;
  using velocity_type = state_type;

  using velocity_policy_t = ExplicitVelocityStandardPolicy<state_type>;

  using type = typename ExplicitImplSelector<tag>::template type<
      true, scalar_type, state_type, system_t, velocity_type, velocity_policy_t
      >;
};

template<class tag, class state_type, class system_t, class policy_t>
struct ExplicitComposeForCustomPolicy
{
  static_assert
  (::pressio::ode::continuous_time_system_with_at_least_velocity<system_t>::value,
   "The system passed to the ExplicitStepper does not meet the required API");

  static_assert
  (::pressio::ode::explicit_state<state_type>::value,
   "Invalid state type for explicit time stepping");

  using scalar_type   = typename ::pressio::traits<state_type>::scalar_type;
  using velocity_type = state_type;
  using time_type = scalar_type;

  static_assert
  (::pressio::ode::explicit_velocity_policy<
      mpl::remove_cvref_t<policy_t>, time_type, state_type, velocity_type, system_t>::value,
   "Invalid rhs policy for explicit time stepping");

  using type = typename ExplicitImplSelector<tag>::template type<
      false, scalar_type, state_type, system_t, velocity_type, policy_t
      >;
};

}}}
#endif  // ODE_EXPLICIT_IMPL_ODE_EXPLICIT_STEPPER_COMPOSE_IMPL_HPP_
