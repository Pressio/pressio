/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_hypred_vel_problem_explicit_step_traits.hpp
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

#ifndef ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_TRAITS_ROM_GALERKIN_HYPRED_VEL_PROBLEM_EXPLICIT_STEP_TRAITS_HPP_
#define ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_TRAITS_ROM_GALERKIN_HYPRED_VEL_PROBLEM_EXPLICIT_STEP_TRAITS_HPP_

namespace pressio{ namespace rom{

//fwd declare problem class
namespace galerkin{ namespace impl{
template <typename ...>
class HypRedVeloProblemExplicitStepContinuousTimeApi;
}}// end namespace pressio::rom::galerkin::impl

namespace details{

template <
  typename stepper_tag,
  typename fom_system_type,
  typename rom_state_type,
  typename decoder_type,
  typename projector_type,
  typename ud_ops_type
  >
struct traits<
  ::pressio::rom::galerkin::impl::HypRedVeloProblemExplicitStepContinuousTimeApi<
    stepper_tag, fom_system_type, rom_state_type, decoder_type, projector_type, ud_ops_type
    >
  >
{
  using common_types_t = ::pressio::rom::galerkin::impl::CommonTraitsContinuousTimeApi<
    stepper_tag, fom_system_type, rom_state_type, decoder_type, ud_ops_type>;

  using fom_system_t		= typename common_types_t::fom_system_t;
  using scalar_t		= typename common_types_t::scalar_t;
  using fom_native_state_t	= typename common_types_t::fom_native_state_t;
  using fom_state_t		= typename common_types_t::fom_state_t;
  using fom_velocity_t		= typename common_types_t::fom_velocity_t;
  using galerkin_state_t	= typename common_types_t::galerkin_state_t;
  using galerkin_native_state_t	= typename common_types_t::galerkin_native_state_t;
  using decoder_t		= typename common_types_t::decoder_t;
  using decoder_jac_t		= typename common_types_t::decoder_jac_t;
  using fom_state_reconstr_t	= typename common_types_t::fom_state_reconstr_t;
  using fom_states_manager_t	= typename common_types_t::fom_states_manager_t;
  using ud_ops_t		= ud_ops_type;
  static constexpr auto binding_sentinel = common_types_t::binding_sentinel;

  // for now, the velocity type is same as state
  using galerkin_velocity_t	= galerkin_state_t;

  using projector_t = projector_type;
  static_assert
  (::pressio::rom::galerkin::constraints::projector_explicit_stepping
   <projector_type, fom_velocity_t, galerkin_velocity_t>::value,
   "Invalid projector passed to Galerkin");

  using velocity_policy_t =
    ::pressio::rom::galerkin::impl::VelocityPolicy<
    galerkin_velocity_t,
    ::pressio::rom::galerkin::impl::Projected<
      projector_t,
      ::pressio::rom::galerkin::impl::FomVelocityPolicy<
	fom_states_manager_t, fom_velocity_t>
      >
    >;

  using stepper_t = ::pressio::ode::ExplicitStepper<
    stepper_tag, galerkin_state_t, fom_system_t, velocity_policy_t>;
};

}}}//end  namespace pressio::rom::galerkin::impl
#endif  // ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_TRAITS_ROM_GALERKIN_HYPRED_VEL_PROBLEM_EXPLICIT_STEP_TRAITS_HPP_
