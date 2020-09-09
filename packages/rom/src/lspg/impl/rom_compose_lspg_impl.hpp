/*
//@HEADER
// ************************************************************************
//
// rom_compose_lspg_impl.hpp
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

#ifndef ROM_LSPG_IMPL_ROM_COMPOSE_LSPG_IMPL_HPP_
#define ROM_LSPG_IMPL_ROM_COMPOSE_LSPG_IMPL_HPP_

#include "./unsteady/continuous_time_api/rom_lspg_unsteady_problem_continuous_time_api.hpp"
#include "./unsteady/discrete_time_api/rom_lspg_unsteady_problem_discrete_time_api.hpp"
#include "./steady/rom_lspg_steady_problem.hpp"

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{

struct Default{};
struct Preconditioned{};
struct Masked{};


//********************
//***** UNSTEADY *****
//********************

template<bool is_probably_discrete_time, typename prob_tag, typename fom_system_t>
struct FindAdapterMistakesUnsteady;

template<typename prob_tag, typename fom_system_t>
struct FindAdapterMistakesUnsteady<true, prob_tag, fom_system_t>
{
  static constexpr bool value =
    ::pressio::rom::what_is_missing_in_passed_discrete_time_system_class<fom_system_t>::value;
};

template<typename fom_system_t>
struct FindAdapterMistakesUnsteady<false, Default, fom_system_t>
{
  static constexpr bool value =
    ::pressio::rom::what_is_missing_in_continuous_time_system_class<fom_system_t>::value;
};

// template<typename fom_system_t>
// struct FindAdapterMistakesUnsteady<false, Preconditioned, fom_system_t>
// {
//   static constexpr bool value =
//     ::pressio::rom::what_is_missing_in_continuous_time_preconditioned_system_class<fom_system_t>::value;
// };


template<typename problem_tag, typename dummy, typename ode_tag, typename fom_system_t, typename ...Args>
struct composeUnsteady
{
  //if we are here, something is wrong, find out what it is

  // check if the stepper tag is wrong
  static_assert
  (::pressio::ode::predicates::is_stepper_tag<ode_tag>::value,
   "\nThe unsteady LSPG problem you are trying to create cannot be created because \
the first template argument is not a valid stepper tag. \
Current choices are: ode::implicitmethods::{Euler, BDF2} for the continuous-time API, \
and  ode::implicitmethods::{Arbitrary} for the discrete-time API.");

  // if we get here, it means that the adapter class is wrong
  // there are two choices: either the user wants to use the discrete-time
  // or the continuous-time API, but we have no way to know here which one
  // they want to use. We also want here to print an error message
  // that is expressive, so to do that we guess what API they want
  // by discriminating based on whether they have a discrete_time_residual_type.
  // We choose this because it is simple enough that we hope the user does not get it wrong.
  // It is very likely that if the system class has that typedef, then
  // the user is trying to use a discrete-time API. If not, then they are
  // using continuous-time API.
  // If the user if actually using the discrete-time API but
  // mispells the "discrete_time_residual_type" typedef, then this breaks.
  // we might need to add logic to handle the case if that typedef is not found
  // because it could mean the user is using the continuous-time api or
  // they mispelled the typedef or they forgot it.

  using error_helper =
    typename std::conditional<
    ::pressio::ode::predicates::has_discrete_time_residual_typedef<fom_system_t>::value,
    FindAdapterMistakesUnsteady<true,  problem_tag, fom_system_t>,
    FindAdapterMistakesUnsteady<false, problem_tag, fom_system_t>
    >::type;
  static_assert(error_helper::value,"");

  // if we get here set to void, but we should never get here in thory
  // because the asserts above should be triggered before
  using type = void;
};


// unsteady default lspg continuous time API
template<typename stepper_tag, typename fom_system_type, typename lspg_state_t, typename ...Args>
struct composeUnsteady<
::pressio::rom::lspg::impl::Default,
mpl::enable_if_t<
::pressio::rom::concepts::continuous_time_implicit_system<fom_system_type>::value
>,
stepper_tag, fom_system_type, lspg_state_t, Args...>
{
  static_assert
  (std::is_same< stepper_tag, ::pressio::ode::implicitmethods::Euler>::value or
   std::is_same< stepper_tag, ::pressio::ode::implicitmethods::BDF2>::value,
   "Unsteady default LSPG with continuous-time API currently only accepts steppers Euler or BDF2");

  using type = ::pressio::rom::lspg::impl::unsteady::ProblemContinuousTimeApi<
            ::pressio::rom::lspg::impl::unsteady::DefaultProblemTraitsContinuousTimeApi,
            stepper_tag, fom_system_type, lspg_state_t, Args...>;
};

// unsteady preconditioned lspg continuous time API
template<typename stepper_tag, typename fom_system_type, typename lspg_state_t, typename ...Args>
struct composeUnsteady<
::pressio::rom::lspg::impl::Preconditioned,
mpl::enable_if_t<
::pressio::rom::concepts::continuous_time_system_preconditionable_rom<fom_system_type>::value
>,
stepper_tag, fom_system_type, lspg_state_t, Args...>
{
  static_assert
  (std::is_same< stepper_tag, ::pressio::ode::implicitmethods::Euler>::value or
   std::is_same< stepper_tag, ::pressio::ode::implicitmethods::BDF2>::value,
   "Unsteady preconditioned LSPG with continuous-time API currently only accepts steppers Euler or BDF2");

  using type = ::pressio::rom::lspg::impl::unsteady::ProblemContinuousTimeApi<
            ::pressio::rom::lspg::impl::unsteady::PreconditionedProblemTraitsContinuousTimeApi,
            stepper_tag, fom_system_type, lspg_state_t, Args...>;
};

// unsteady masked lspg continuous time API
template<typename stepper_tag, typename fom_system_type, typename lspg_state_t, typename ...Args>
struct composeUnsteady<
::pressio::rom::lspg::impl::Masked,
mpl::enable_if_t<
::pressio::rom::concepts::continuous_time_system_maskable_rom<fom_system_type>::value
>,
stepper_tag, fom_system_type, lspg_state_t, Args...>
{
  static_assert
  (std::is_same< stepper_tag, ::pressio::ode::implicitmethods::Euler>::value or
   std::is_same< stepper_tag, ::pressio::ode::implicitmethods::BDF2>::value,
   "Unsteady masked LSPG with continuous-time API currently only accepts steppers Euler or BDF2");

  using type = ::pressio::rom::lspg::impl::unsteady::ProblemContinuousTimeApi<
            ::pressio::rom::lspg::impl::unsteady::MaskedProblemTraitsContinuousTimeApi,
            stepper_tag, fom_system_type, lspg_state_t, Args...>;
};


// unsteady default lspg discrete time api
template<typename stepper_tag, typename fom_system_type, typename lspg_state_t, typename ...Args>
struct composeUnsteady<
::pressio::rom::lspg::impl::Default,
mpl::enable_if_t<
::pressio::rom::concepts::discrete_time_system<fom_system_type>::value
>,
stepper_tag, fom_system_type, lspg_state_t, Args...>
{
  static_assert
  (std::is_same< stepper_tag, ::pressio::ode::implicitmethods::Arbitrary>::value,
   "Unsteady default LSPG with discrete-time API currently only accepts Arbitrary stepper");

  using type = ::pressio::rom::lspg::impl::unsteady::ProblemDiscreteTimeApi<
            ::pressio::rom::lspg::impl::unsteady::DefaultProblemTraitsDiscreteTimeApi,
            stepper_tag, fom_system_type, lspg_state_t, Args...>;
};



//********************
//***** STEADY *****
//********************

template<typename problem_tag, typename fom_system_t, typename ...Args>
struct composeSteady
{
  //if we are here, something is wrong, find out if it is the
  // API of the system_type or something else

//   static_assert
//   (::pressio::rom::concepts::continuous_time_implicit_system<fom_system_t>::value
//    or ::pressio::rom::concepts::discrete_time_system<fom_system_t>::value,
//    "\nThe LSPG problem you are trying to create cannot be created because \
// the fom adapter class you are using does not meet neither the \
// continuous-time nor the discrete-time API. \n\
// If you are tring to use the discrete-time API, \
// you can call right after your declare your fom adapter type \
// the following function which will generate a compile-time message: \
// constexpr auto dummy = pressio::rom::what_is_missing_in_my_discrete_time_system_class<your_adapter_type>();");

  using type = void;
};


// default
template<typename fom_system_type, typename lspg_state_t, typename ...Args>
struct composeSteady<
::pressio::rom::lspg::impl::Default,
mpl::enable_if_t<
::pressio::rom::concepts::steady_system<fom_system_type>::value
>,
fom_system_type, lspg_state_t, Args...>
{
  using type = ::pressio::rom::lspg::impl::steady::ProblemSteady<
            ::pressio::rom::lspg::impl::steady::DefaultProblemTraits,
            fom_system_type, lspg_state_t, Args...>;
};

// preconditionale rom
template<typename fom_system_type, typename lspg_state_t, typename ...Args>
struct composeSteady<
::pressio::rom::lspg::impl::Preconditioned,
mpl::enable_if_t<
::pressio::rom::concepts::steady_system_preconditionable_rom<fom_system_type>::value
>,
fom_system_type, lspg_state_t, Args...>
{
  using type = ::pressio::rom::lspg::impl::steady::ProblemSteady<
            ::pressio::rom::lspg::impl::steady::PreconditionedProblemTraits,
            fom_system_type, lspg_state_t, Args...>;
};

}}}}
#endif  // ROM_LSPG_IMPL_ROM_COMPOSE_LSPG_IMPL_HPP_
