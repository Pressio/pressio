/*
//@HEADER
// ************************************************************************
//
// rom_lspg_type_generator_default.hpp
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

#ifndef ROM_LSPG_STEADY_TYPE_GENERATOR_DEFAULT_HPP_
#define ROM_LSPG_STEADY_TYPE_GENERATOR_DEFAULT_HPP_

#include "../../rom_lspg_type_generator_common.hpp"

namespace pressio{ namespace rom{

template <
  typename fom_type,
  typename decoder_type,
  typename lspg_state_type
  >
struct DefaultLSPGSteadyTypeGenerator
  : LSPGCommonTypes<
  fom_type, decoder_type, lspg_state_type
  >{

  using this_t = DefaultLSPGSteadyTypeGenerator
    <fom_type, decoder_type, lspg_state_type>;

  using base_t = LSPGCommonTypes
    <fom_type, decoder_type, lspg_state_type>;

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

  static constexpr bool is_steady = true;

  /* lspg_matrix_t is type of J*decoder_jac_t (in the most basic case) where
   * * J is the jacobian of the fom rhs
   * * decoder_jac_t is the type of the decoder jacobian
   * In more complex cases, we might have (something)*J*decoder_jac_t,
   * where (something) is product of few matrices.
   * For now, set lspg_matrix_t to be of same type as decoder_jac_t
   * if phi is MV<>, then lspg_matrix_t = containers::MV<>
   * if phi is Matrix<>, then we have containers::Matrix<>
   * not a bad assumption since all matrices are left-applied to decoder_jac_t
   */
  using lspg_matrix_t		= decoder_jac_t;

  // policy for evaluating the rhs of the fom object (<true> for steady overload)
  using fom_eval_rhs_policy_t	= ::pressio::rom::policy::EvaluateFomVelocityDefault<this_t::is_steady>;

  // policy for left multiplying the fom jacobian with decoder_jac_t
  // possibly involving other stuff like explained above (<true> for steady overload
  using fom_apply_jac_policy_t	= ::pressio::rom::policy::ApplyFomJacobianDefault<this_t::is_steady>;

  // Policy defining how to compute the LSPG residual
  using lspg_residual_policy_t	= ::pressio::rom::LSPGSteadyResidualPolicy<
	fom_states_data, fom_velocity_data, fom_eval_rhs_policy_t>;

  // policy defining how to compute the LSPG jacobian
  using lspg_jacobian_policy_t	= ::pressio::rom::LSPGSteadyJacobianPolicy<
    fom_states_data, lspg_matrix_t, fom_apply_jac_policy_t, decoder_t>;

  // system's type
  using lspg_system_t		= ::pressio::rom::LSPGSteadySystem<
    fom_t, lspg_state_type, lspg_residual_t, lspg_matrix_t,
    lspg_residual_policy_t, lspg_jacobian_policy_t>;

};//end class

}}//end  namespace pressio::rom
#endif
