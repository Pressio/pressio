/*
//@HEADER
// ************************************************************************
//
// ode_explicit_stepper_impl.hpp
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

#ifndef ODE_EXPLICIT_IMPL_ODE_EXPLICIT_STEPPER_IMPL_HPP_
#define ODE_EXPLICIT_IMPL_ODE_EXPLICIT_STEPPER_IMPL_HPP_

#include "ode_explicit_euler_stepper_impl.hpp"
#include "ode_explicit_runge_kutta4_stepper_impl.hpp"

namespace pressio{ namespace ode{ namespace explicitmethods{ namespace impl{

template<typename tag, typename ...Args>
struct UserDefinedOpsFilter{ using type = void; };

template<typename scalar_t, typename state_t, typename velocity_t, typename ...Args>
struct UserDefinedOpsFilter<
::pressio::ode::explicitmethods::Euler, scalar_t, state_t, velocity_t, Args...
>
{
  // check if user passed an ops
  using ic4 = ::pressio::mpl::variadic::find_if_quaternary_pred_t<
    scalar_t, state_t, velocity_t,
    ::pressio::ode::concepts::user_defined_ops_for_explicit_euler, Args...>;
  using type = ::pressio::mpl::variadic::at_or_t<void, ic4::value, Args...>;
};

template<typename scalar_t, typename state_t, typename velocity_t, typename ...Args>
struct UserDefinedOpsFilter<
::pressio::ode::explicitmethods::RungeKutta4, scalar_t, state_t, velocity_t, Args...
>
{
  // check if user passed an ops
  using ic4 = ::pressio::mpl::variadic::find_if_quaternary_pred_t<
    scalar_t, state_t, velocity_t,
    ::pressio::ode::concepts::user_defined_ops_for_explicit_rk4, Args...>;
  using type = ::pressio::mpl::variadic::at_or_t<void, ic4::value, Args...>;
};


template<typename tag>
struct ImplSelector{
	template<typename ...Args> using type = void;
};
template<>
struct ImplSelector<::pressio::ode::explicitmethods::Euler>{
	template<typename ...Args> 
	using type = ::pressio::ode::explicitmethods::impl::ExplicitEulerStepperImpl<Args...>;
};
template<>
struct ImplSelector<::pressio::ode::explicitmethods::RungeKutta4>{
  template<typename ...Args> 
  using type = ::pressio::ode::explicitmethods::impl::ExplicitRungeKutta4StepperImpl<Args...>;
};


// Stepper< tag, state_type, system, Args...>
// Args can contain : ud_ops, velocity_policy, velocity_type in any order

template<typename tag, typename state_type, typename system_t, typename ...Args>
struct Stepper
{ 
  using state_t	= state_type;

  // // check if scalar is provided in Args
  // using ic0 = ::pressio::mpl::variadic::find_if_unary_pred_t<
  //   std::is_floating_point, Args...>;
  // using scalar_t = ::pressio::mpl::variadic::at_or_t<
  //   void, ic0::value, Args...>;
  // static_assert(std::is_void<scalar_t>::value == false,
		//  "You need a scalar_type in the ExplicitStepper templates");
  using scalar_t = typename ::pressio::containers::details::traits<state_type>::scalar_t;

  static_assert(::pressio::ode::concepts::continuous_time_explicit_system<system_t>::value,
       "Invalid system passed to the ExplicitStepper");

  // check args for a valid velocity type
  using ic2 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::ode::concepts::explicit_velocity, Args...>;
  // if a velocity type is NOT found, then set it equal to the state
  using velocity_t = ::pressio::mpl::variadic::at_or_t<state_t, ic2::value, Args...>;
  static_assert(std::is_void<system_t>::value == false,
  		 "The velicity type cannot be void");

  // this is the standard velocity policy (just typedef, it is only used
  // if the user does not pass a user-defined policy)
  using standard_velocity_policy_t = ::pressio::ode::explicitmethods::policy::VelocityStandardPolicy<
    state_t, system_t, velocity_t>;

  // check Args if a user-defined velocity policy is passed
  using ic3 = ::pressio::mpl::variadic::find_if_quinary_pred_t<
    scalar_t, state_type, velocity_t, system_t,
    ::pressio::ode::concepts::explicit_velocity_policy, Args...>;
  using velocity_policy_t = ::pressio::mpl::variadic::at_or_t<standard_velocity_policy_t, ic3::value, Args...>;

  // check for user-defined ops
  using ops_t = typename UserDefinedOpsFilter<tag, scalar_t, state_t, velocity_t, Args...>::type;

  // implementation class type
  using type = typename ImplSelector<tag>::template type<scalar_t, state_t, system_t, velocity_t, 
  			velocity_policy_t, standard_velocity_policy_t, ops_t>;
};

}}}} // end namespace pressio::ode::explicitmethods::impl
#endif  // ODE_EXPLICIT_IMPL_ODE_EXPLICIT_STEPPER_IMPL_HPP_
