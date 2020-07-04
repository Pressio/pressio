/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_problem_type_generator_default_velocity_api.hpp
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

#ifndef ROM_GALERKIN_PROBLEM_TYPE_GENERATOR_DEFAULT_VELOCITY_API_HPP_
#define ROM_GALERKIN_PROBLEM_TYPE_GENERATOR_DEFAULT_VELOCITY_API_HPP_

#include "rom_galerkin_type_generator_common_velocity_api.hpp"
#include "rom_galerkin_explicit_velocity_policy.hpp"

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template <
  typename fom_type,
  typename stepper_tag,
  typename rom_state_type,
  typename ...Args
  >
struct DefaultProblemTypeGeneratorVelocityApi
{
  /* here, the fom_type must satisfy the velocity api */
  static_assert( ::pressio::rom::meta::model_meets_velocity_api_for_galerkin<fom_type>::value,
		 "\nUsing DefaultProblemTypeGeneratorVelocityApi \n \
requires a fom adapter class that meets the velocity api. \n \
However, the fom/adapter type you passed does not meet the velocity api. \n \
Verify the fom/adapter class you are using meets the velocity api.");

  static_assert( std::is_same< stepper_tag,  ::pressio::ode::explicitmethods::Euler>::value or
		 std::is_same< stepper_tag,  ::pressio::ode::explicitmethods::RungeKutta4>::value,
		 "Velocity-API-based Galerkin currently only supports explicit time-stepping. \n \
You need to pass a valid stepper_tag from ::pressio::ode::explicitmethods");

  // pick the common types holder
  using common_types_t = GalerkinCommonTypesVelocityApi<fom_type, rom_state_type, Args...>;

  using fom_t			= typename common_types_t::fom_t;
  using scalar_t		= typename common_types_t::scalar_t;
  using fom_native_state_t	= typename common_types_t::fom_native_state_t;
  using fom_state_t		= typename common_types_t::fom_state_t;
  using fom_velocity_t		= typename common_types_t::fom_velocity_t;
  using galerkin_state_t	= typename common_types_t::galerkin_state_t;
  using galerkin_native_state_t	= typename common_types_t::galerkin_native_state_t;
  using galerkin_residual_t	= typename common_types_t::galerkin_residual_t;
  using decoder_t		= typename common_types_t::decoder_t;
  using decoder_jac_t		= typename common_types_t::decoder_jac_t;
  using fom_state_reconstr_t	= typename common_types_t::fom_state_reconstr_t;
  using fom_states_manager_t		= typename common_types_t::fom_states_manager_t;
  using ud_ops_t		= typename common_types_t::ud_ops_t;

  // policy for evaluating the ode velocity
  using velocity_policy_t =
    ::pressio::rom::galerkin::impl::ExplicitVelocityPolicy<rom_state_type, 
    fom_states_manager_t, fom_velocity_t, decoder_t, ud_ops_t>;

  // declare type of stepper object
  using stepper_t = ::pressio::ode::explicitmethods::Stepper<
    stepper_tag, rom_state_type, fom_t, galerkin_residual_t, velocity_policy_t, scalar_t>;

};//end class

}}}}//end  namespace pressio::rom::galerkin::impl
#endif
