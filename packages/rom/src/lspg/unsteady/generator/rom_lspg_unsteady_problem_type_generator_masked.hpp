/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_problem_type_generator_masked.hpp
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

#ifndef ROM_LSPG_TYPE_GENERATOR_MASKED_HPP_
#define ROM_LSPG_TYPE_GENERATOR_MASKED_HPP_

#include "../../rom_lspg_type_generator_common.hpp"

namespace pressio{ namespace rom{

template <
  typename fom_type,
  ode::ImplicitEnum odeName,
  typename decoder_type,
  typename lspg_state_type,
  typename ud_ops = void
  >
struct MaskedLSPGTypeGenerator
  : LSPGCommonTypes<
  fom_type, decoder_type, lspg_state_type, odeName, ud_ops
  >{

  using base_t = LSPGCommonTypes<
    fom_type, decoder_type, lspg_state_type, odeName, ud_ops>;

  using typename base_t::fom_t;
  using typename base_t::scalar_t;
  using typename base_t::fom_native_state_t;
  using typename base_t::fom_state_t;
  using typename base_t::fom_velocity_t;
  using typename base_t::lspg_state_t;
  using typename base_t::lspg_residual_t;
  using typename base_t::decoder_t;
  using typename base_t::decoder_jac_t;
  using typename base_t::fom_state_reconstr_t;
  using typename base_t::fom_states_data;
  using typename base_t::fom_velocity_data;
  using typename base_t::ud_ops_t;

  /* lspg_matrix_t is type of J*decoder_jac_t (in the most basic case) where
   * * J is the jacobian of the fom rhs
   * * decoder_jac_t is the type of the decoder jacobian
   * In more complex cases, we might have (something)*J*decoder_jac_t,
   * where (something) is product of few matrices.
   * For now, set lspg_matrix_t to be of same type as decoder_jac_t
   * if phi is MV<>, then lspg_matrix_t = containers::MV<>
   * if phi is Matrix<>, then we have containers::Matrix<>
   * not bad assumption since all matrices are left-applied to decoder_jac_t
   */
  using lspg_matrix_t		= decoder_jac_t;

  // policy for evaluating the rhs of the fom object (<false> for unsteady overload)
  using fom_eval_velocity_policy_t	= ::pressio::rom::policy::EvaluateFomVelocityDefault<false>;

  // policy for left multiplying the fom jacobian with decoder_jac_t
  // possibly involving other stuff like explained above (<false> for unsteady overload)
  using fom_apply_jac_policy_t	= ::pressio::rom::policy::ApplyFomJacobianDefault<false>;

  // policy defining how to compute the LSPG time-discrete residual
  using lspg_residual_policy_t =
    rom::decorator::Masked<
    rom::LSPGResidualPolicy<
      fom_states_data, fom_velocity_data, fom_eval_velocity_policy_t, ud_ops
      >
    >;

  // policy defining how to compute the LSPG time-discrete jacobian
  using lspg_jacobian_policy_t	=
    rom::decorator::Masked<
    rom::LSPGJacobianPolicy<
      fom_states_data, lspg_matrix_t, fom_apply_jac_policy_t, decoder_t, ud_ops
      >
    >;

  // auxiliary stepper
  using aux_stepper_t = typename auxStepperHelper<
    odeName, lspg_state_type,
    lspg_residual_t, lspg_matrix_t,
    fom_type, lspg_residual_policy_t,
    lspg_jacobian_policy_t, scalar_t>::type;

  // stepper object type
  using lspg_stepper_t		= ode::ImplicitStepper<
    odeName, lspg_state_type,
    lspg_residual_t, lspg_matrix_t,
    fom_type, aux_stepper_t,
    lspg_residual_policy_t, lspg_jacobian_policy_t, scalar_t>;

};//end class


}}//end  namespace pressio::rom
#endif
