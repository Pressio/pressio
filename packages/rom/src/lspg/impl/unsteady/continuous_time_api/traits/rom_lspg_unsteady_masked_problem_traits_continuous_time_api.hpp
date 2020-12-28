/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_masked_problem_traits_continuous_time_api.hpp
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

#ifndef ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_TRAITS_ROM_LSPG_UNSTEADY_MASKED_PROBLEM_TRAITS_CONTINUOUS_TIME_API_HPP_
#define ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_TRAITS_ROM_LSPG_UNSTEADY_MASKED_PROBLEM_TRAITS_CONTINUOUS_TIME_API_HPP_

namespace pressio{ namespace rom{

namespace lspg{ namespace impl{ namespace unsteady{
template <typename...>
class MaskedProblemContinuousTimeApi;
}}}// end namespace pressio::rom::lspg::impl::unsteady

namespace details{

template <
  typename stepper_tag,
  typename fom_system_type,
  typename lspg_state_type,
  typename decoder_type,
  typename masker_type,
  typename ud_ops_type
  >
struct traits<
  ::pressio::rom::lspg::impl::unsteady::MaskedProblemContinuousTimeApi<
    stepper_tag, fom_system_type, lspg_state_type, decoder_type, masker_type, ud_ops_type
    >
  >
{
  static const bool is_steady_lspg = false;
  static const bool is_unsteady_lspg = true;

  using common_types_t =
    ::pressio::rom::lspg::impl::unsteady::CommonTraitsContinuousTimeApi<
    stepper_tag, fom_system_type, lspg_state_type, decoder_type, ud_ops_type>;

  using fom_system_t		= typename common_types_t::fom_system_t;
  using scalar_t		= typename common_types_t::scalar_t;
  using fom_native_state_t	= typename common_types_t::fom_native_state_t;
  using fom_state_t		= typename common_types_t::fom_state_t;
  using fom_velocity_t		= typename common_types_t::fom_velocity_t;
  using lspg_state_t		= typename common_types_t::lspg_state_t;
  using lspg_native_state_t	= typename common_types_t::lspg_native_state_t;
  using lspg_residual_t		= typename common_types_t::lspg_residual_t;
  using decoder_t		= typename common_types_t::decoder_t;
  using decoder_jac_t		= typename common_types_t::decoder_jac_t;
  using lspg_jacobian_t		= typename common_types_t::lspg_jacobian_t;
  using fom_state_reconstr_t	= typename common_types_t::fom_state_reconstr_t;
  using fom_states_manager_t	= typename common_types_t::fom_states_manager_t;
  using ud_ops_t = ud_ops_type;
  static constexpr auto binding_sentinel = common_types_t::binding_sentinel;

  using masker_t = masker_type;
  static_assert
  (
   ::pressio::rom::lspg::constraints::unsteady_masker<
   masker_t,
   scalar_t,
   typename ::pressio::containers::details::traits<fom_velocity_t>::wrapped_t,
   typename ::pressio::containers::details::traits<lspg_jacobian_t>::wrapped_t
   >::value,
   "Invalid masker type passed to unsteady LSPG"
  );

  using residual_policy_t =
    ::pressio::rom::lspg::decorator::Masked<
    masker_t,
    ::pressio::rom::lspg::impl::unsteady::ResidualPolicyContinuousTimeApi<
      lspg_residual_t, fom_states_manager_t, ud_ops_t
      >
    >;

  using jacobian_policy_t	=
    ::pressio::rom::lspg::decorator::Masked<
    masker_t,
    ::pressio::rom::lspg::impl::unsteady::JacobianPolicyContinuousTimeApi<
      fom_states_manager_t, lspg_jacobian_t, decoder_t, ud_ops_t
      >
    >;

  // auxiliary stepper
  using aux_stepper_t = typename ::pressio::rom::impl::auxiliaryStepperHelper<
    stepper_tag, lspg_state_type, lspg_residual_t, lspg_jacobian_t, fom_system_type,
    residual_policy_t, jacobian_policy_t>::type;

  // stepper object type
  using stepper_t = ::pressio::ode::ImplicitStepper<
    stepper_tag, lspg_state_type, lspg_residual_t, lspg_jacobian_t, fom_system_type,
    aux_stepper_t, residual_policy_t, jacobian_policy_t>;

};//end class

}}}//end  namespace pressio::rom::lspg::unstedy::impl
#endif  // ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_TRAITS_ROM_LSPG_UNSTEADY_MASKED_PROBLEM_TRAITS_CONTINUOUS_TIME_API_HPP_
