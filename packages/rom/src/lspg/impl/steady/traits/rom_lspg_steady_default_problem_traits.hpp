/*
//@HEADER
// ************************************************************************
//
// rom_lspg_steady_default_problem_traits.hpp
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

#ifndef ROM_LSPG_IMPL_STEADY_TRAITS_ROM_LSPG_STEADY_DEFAULT_PROBLEM_TRAITS_HPP_
#define ROM_LSPG_IMPL_STEADY_TRAITS_ROM_LSPG_STEADY_DEFAULT_PROBLEM_TRAITS_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{ namespace steady{

template <
  typename fom_type,
  typename lspg_state_type,
  typename decoder_type,
  typename ...Args
  >
struct DefaultProblemTraits
  : CommonTraits<fom_type, decoder_type, lspg_state_type, Args...>
{

  using base_t = ::pressio::rom::lspg::impl::steady::CommonTraits<fom_type, decoder_type, lspg_state_type, Args...>;

  using typename base_t::fom_t;
  using typename base_t::scalar_t;
  using typename base_t::fom_native_state_t;
  using typename base_t::fom_state_t;
  using typename base_t::fom_residual_t;
  using typename base_t::lspg_state_t;
  using typename base_t::lspg_residual_t;
  using typename base_t::decoder_t;
  using typename base_t::decoder_jac_t;
  using typename base_t::fom_state_reconstr_t;
  using typename base_t::fom_states_manager_t;
  using typename base_t::ud_ops_t;

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

  // Policy defining how to compute the LSPG residual
  using lspg_residual_policy_t	= ::pressio::rom::lspg::impl::steady::ResidualPolicy<
	lspg_residual_t, fom_states_manager_t, void>;

  // policy defining how to compute the LSPG jacobian
  using lspg_jacobian_policy_t	= ::pressio::rom::lspg::impl::steady::JacobianPolicy<
    fom_states_manager_t, lspg_matrix_t, decoder_t, void>;

  // system's type
  using lspg_system_t		= ::pressio::rom::lspg::impl::steady::System<
    fom_t, lspg_state_type, lspg_residual_t, lspg_matrix_t,
    lspg_residual_policy_t, lspg_jacobian_policy_t>;

};//end class

}}}}}
#endif  // ROM_LSPG_IMPL_STEADY_TRAITS_ROM_LSPG_STEADY_DEFAULT_PROBLEM_TRAITS_HPP_
