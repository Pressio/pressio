/*
//@HEADER
// ************************************************************************
//
// rom_compose_galerkin_impl.hpp
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

#ifndef ROM_GALERKIN_IMPL_ROM_COMPOSE_GALERKIN_IMPL_HPP_
#define ROM_GALERKIN_IMPL_ROM_COMPOSE_GALERKIN_IMPL_HPP_

// continuous time API
#include "./continuous_time_api/rom_galerkin_explicit_velocity_policy.hpp"
#include "./continuous_time_api/traits/rom_galerkin_common_traits_continuous_time_api.hpp"
#include "./continuous_time_api/traits/rom_galerkin_default_problem_traits_continuous_time_api.hpp"
#include "./continuous_time_api/rom_galerkin_default_problem_continuous_time_api.hpp"

// discrete time API
#include "./discrete_time_api/policies/rom_galerkin_residual_policy_discrete_time_api.hpp"
#include "./discrete_time_api/policies/rom_galerkin_jacobian_policy_discrete_time_api.hpp"
#include "./discrete_time_api/traits/rom_galerkin_common_traits_discrete_time_api.hpp"
#include "./discrete_time_api/traits/rom_galerkin_default_problem_traits_discrete_time_api.hpp"
#include "./discrete_time_api/rom_galerkin_default_problem_discrete_time_api.hpp"

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

// tags for the Galerkin problem
struct Default{};

//---------------------------------------------
// helpers classes to print
// errors in the fom adapter at compiletime
//---------------------------------------------
template<bool is_probably_discrete_time, typename prob_tag, typename fom_system_t>
struct FindAdapterMistakes;

template<typename prob_tag, typename fom_system_t>
struct FindAdapterMistakes<true, prob_tag, fom_system_t>
{
  static constexpr bool value =
    ::pressio::rom::find_discrepancies_with_discrete_time_system_api<fom_system_t>::value;
};

template<typename fom_system_t>
struct FindAdapterMistakes<false, Default, fom_system_t>
{
  static constexpr bool value =
    ::pressio::rom::find_discrepancies_with_continuous_time_explicit_system_api<fom_system_t>::value;
};

// //------------------------
// // specialize compose
// //------------------------
template<
  typename problem_tag,
  typename dummy,
  typename ode_tag,
  typename fom_system_t,
  typename ...Args>
struct compose
{
  //if we get here, something is wrong, find out what it is

  // check if the stepper tag is wrong
  static_assert
  (::pressio::ode::predicates::is_stepper_tag<ode_tag>::value,
   "\nThe Galerkin problem you are trying to create cannot be created because \
the first template argument is not a valid stepper tag.			\
Valid choices: ode::explicitmethods::{Euler, RungeKutta4} for the continuous-time API, \
and ode::implicitmethods::{Arbitrary} for the discrete-time API.");

/*  if here, it means that the adapter class is the problem.
    There are two scenarios: either the user wants to use the discrete-time
    or the continuous-time API. But we have no way to know here which one
    they WANT to use and want to print error messages detailed enough.
    So to do that we guess what API they want
    by discriminating on whether the adapter has a discrete_time_residual_type.
    We choose this because it is simple enough that we hope the user does not get it wrong.
    We believe that if the system class has the discrete_time_residual_type typedef,
    then it is very likely the user is trying to use a discrete-time API.
    If not, then they are trying to use the continuous-time API.
    If the user if actually using the discrete-time API but
    mispells or forgets the "discrete_time_residual_type" typedef,
    then the logic here breaks, and we might need to add logic to handle the
    case if that typedef is not found, because it could mean the user is
    using the continuous-time api or mispelled the typedef or forgot it. */

  using helper =
    typename std::conditional<
    ::pressio::ode::predicates::has_discrete_time_residual_typedef<fom_system_t>::value,
    FindAdapterMistakes<true , problem_tag, fom_system_t>,
    FindAdapterMistakes<false, problem_tag, fom_system_t>
    >::type;
  // I need this static assert otherwise it won't trigger errors
  static_assert(helper::value,"");

  // we should never get here, things should always fail before
  using type = void;
};

///////////////////////////////////////////////////////
//////       CONTINUOUS TIME API
///////////////////////////////////////////////////////
/**** 
  default Galerkin, pressio ops
***/
template<
  typename stepper_tag, 
  typename fom_system_type, 
  typename galerkin_state_type, 
  typename decoder_type
  >
struct compose<
::pressio::rom::galerkin::impl::Default,
mpl::enable_if_t<
::pressio::rom::concepts::continuous_time_explicit_system<fom_system_type>::value
>,
stepper_tag, fom_system_type, galerkin_state_type, decoder_type>
{
  static_assert
  (std::is_same< stepper_tag, ::pressio::ode::explicitmethods::Euler>::value or
   std::is_same< stepper_tag, ::pressio::ode::explicitmethods::RungeKutta4>::value,
   "Default Galerkin with continuous-time API currently only accepts Forward Euler or RK4");

  using type = ::pressio::rom::galerkin::impl::DefaultProblemContinuousTimeApi<
    stepper_tag, fom_system_type, galerkin_state_type, decoder_type, void>;
};

/**** 
  default Galerkin, user-defined ops
***/
template<
  typename stepper_tag, 
  typename fom_system_type, 
  typename galerkin_state_type, 
  typename decoder_type,
  typename ud_ops_type
  >
struct compose<
::pressio::rom::galerkin::impl::Default,
mpl::enable_if_t<
::pressio::rom::concepts::continuous_time_explicit_system<fom_system_type>::value
>,
stepper_tag, fom_system_type, galerkin_state_type, decoder_type, ud_ops_type>
{
  static_assert
  (std::is_same< stepper_tag, ::pressio::ode::explicitmethods::Euler>::value or
   std::is_same< stepper_tag, ::pressio::ode::explicitmethods::RungeKutta4>::value,
   "Default Galerkin with continuous-time API currently only accepts Forward Euler or RK4");

  using type = ::pressio::rom::galerkin::impl::DefaultProblemContinuousTimeApi<
    stepper_tag, fom_system_type, galerkin_state_type, decoder_type, ud_ops_type>;
};


///////////////////////////////////////////////////////
//////       DISCRETE TIME API
///////////////////////////////////////////////////////
/**** 
  default Galerkin, pressio ops
***/
template<
  typename stepper_tag, 
  typename fom_system_type, 
  typename galerkin_state_type, 
  typename rom_jacobian_type,
  typename decoder_type,
  typename ...Args
  >
struct compose<
::pressio::rom::galerkin::impl::Default,
mpl::enable_if_t<
::pressio::rom::concepts::discrete_time_system<fom_system_type>::value
>,
stepper_tag, fom_system_type, galerkin_state_type, rom_jacobian_type, decoder_type, Args...>
{
  static_assert
  (std::is_same< stepper_tag, ::pressio::ode::implicitmethods::Arbitrary>::value,
   "Default Galerkin with discrete-time API currently only accepts Arbitrary stepper");

  using type = ::pressio::rom::galerkin::impl::DefaultProblemDiscreteTimeApi<
      stepper_tag, fom_system_type, galerkin_state_type,  
      rom_jacobian_type, decoder_type, Args...>;
};

}}}}
#endif  // ROM_GALERKIN_IMPL_ROM_COMPOSE_GALERKIN_IMPL_HPP_
