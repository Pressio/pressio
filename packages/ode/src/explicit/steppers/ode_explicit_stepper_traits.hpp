/*
//@HEADER
// ************************************************************************
//
// ode_explicit_stepper_traits.hpp
//                     		      Pressio 
// Copyright 2019 National Technology & Engineering Solutions of Sandia,LLC 
//							      (NTESS)
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

#ifndef ODE_STEPPERS_EXPLICIT_STEPPERS_EXPLICIT_STEPPER_TRAITS_HPP_
#define ODE_STEPPERS_EXPLICIT_STEPPERS_EXPLICIT_STEPPER_TRAITS_HPP_

#include "../../ode_fwd.hpp"

namespace pressio{ namespace ode{ namespace details{

/*
 * Euler
 */
template<
  typename state_type, typename ...Args
  >
struct traits<
  ExplicitStepper<ExplicitEnum::Euler, state_type, Args...>
  >{

  static constexpr bool is_implicit = false;
  static constexpr bool is_explicit = true;
  using order_t = unsigned int;
  static constexpr order_t order_value = 1;

  using state_t	   = state_type;

  // check if scalar is provided in Args
  using ic0 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    std::is_floating_point, Args...>;
  using scalar_t = ::pressio::mpl::variadic::at_or_t<
    void, ic0::value, Args...>;
  static_assert( std::is_void<scalar_t>::value == false,
		 "You need a scalar_type in the ExplicitStepper templates");

  // check args for a valid model type
  using ic1 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::ode::meta::is_legitimate_model_for_explicit_ode, Args...>;
  using model_t = ::pressio::mpl::variadic::at_or_t<void, ic1::value, Args...>;
  static_assert( std::is_void<model_t>::value == false,
  		 "Non-admissible model type argument passed to the ExplicitStepper");

  // check args for a valid velocity type
  using ic2 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::ode::meta::is_legitimate_explicit_velocity_type, Args...>;
  // if a velocity type is NOT found, then set it equal to the state
  using velocity_t = ::pressio::mpl::variadic::at_or_t<state_t, ic2::value, Args...>;
  static_assert( std::is_void<model_t>::value == false,
  		 "The velicity type cannot be void");

  // this is the standard velocity policy (just typedef, it is only used
  // if the user does not pass a user-defined policy)
  using standard_res_policy_t = policy::ExplicitVelocityStandardPolicy<
    state_t, model_t, velocity_t>;

  // check Args if a user-defined velocity policy is passed
  using ic3 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::ode::meta::is_legitimate_explicit_velocity_policy, Args...>;
  using velocity_policy_t = ::pressio::mpl::variadic::at_or_t
    <standard_res_policy_t, ic3::value, Args...>;

  // check if user passed an ops
  using ic4 = ::pressio::mpl::variadic::find_if_quaternary_pred_t<
    scalar_t, state_t, velocity_t,
    ::pressio::ode::meta::is_valid_user_defined_ops_for_explicit_euler, Args...>;
  using ops_t = ::pressio::mpl::variadic::at_or_t<void, ic4::value, Args...>;

  // implementation class type
  using impl_t = impl::ExplicitEulerStepperImpl
    <scalar_t, state_t, model_t, velocity_t, velocity_policy_t, ops_t>;
};


/*
 * RK4
 */
template<
  typename state_type, typename ...Args
  >
struct traits<
  ExplicitStepper<ExplicitEnum::RungeKutta4, state_type, Args...>
  >{

  static constexpr bool is_implicit = false;
  static constexpr bool is_explicit = true;
  using order_t = unsigned int;
  static constexpr order_t order_value = 4;

  using state_t	   = state_type;

  // check if scalar is provided in Args
  using ic0 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    std::is_floating_point, Args...>;
  using scalar_t = ::pressio::mpl::variadic::at_or_t<
    void, ic0::value, Args...>;
  static_assert( std::is_void<scalar_t>::value == false,
		 "You need a scalar_type in the ExplicitStepper templates");

  // check args for a valid model type
  using ic1 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::ode::meta::is_legitimate_model_for_explicit_ode, Args...>;
  using model_t = ::pressio::mpl::variadic::at_or_t<void, ic1::value, Args...>;
  static_assert( std::is_void<model_t>::value == false,
  		 "Non-admissible model type argument passed to the ExplicitStepper");

  // check args for a valid velocity type
  using ic2 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::ode::meta::is_legitimate_explicit_velocity_type, Args...>;
  // if a velocity type is NOT found, then set it equal to the state
  using velocity_t = ::pressio::mpl::variadic::at_or_t<state_t, ic2::value, Args...>;
  static_assert( std::is_void<model_t>::value == false,
  		 "The velicity type cannot be void");

  // this is the standard velocity policy (just typedef, it is only used
  // if the user does not pass a user-defined policy)
  using standard_res_policy_t = policy::ExplicitVelocityStandardPolicy<
    state_t, model_t, velocity_t>;

  // check Args if a user-defined velocity policy is passed
  using ic3 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::ode::meta::is_legitimate_explicit_velocity_policy, Args...>;
  using velocity_policy_t = ::pressio::mpl::variadic::at_or_t
    <standard_res_policy_t, ic3::value, Args...>;

  // check if user passed an ops
  using ic4 = ::pressio::mpl::variadic::find_if_quaternary_pred_t<
    scalar_t, state_t, velocity_t,
    ::pressio::ode::meta::is_valid_user_defined_ops_for_explicit_rk4, Args...>;
  using ops_t = ::pressio::mpl::variadic::at_or_t<void, ic4::value, Args...>;

  // implementation class type
  using impl_t = impl::ExplicitRungeKutta4StepperImpl
    <scalar_t, state_t, model_t, velocity_t, velocity_policy_t, ops_t>;
};


}}}//end namespace pressio::ode::details
#endif
