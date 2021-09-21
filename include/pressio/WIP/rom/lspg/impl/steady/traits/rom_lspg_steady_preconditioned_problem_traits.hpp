/*
//@HEADER
// ************************************************************************
//
// rom_lspg_steady_preconditioned_problem_traits.hpp
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

#ifndef ROM_LSPG_IMPL_STEADY_TRAITS_ROM_LSPG_STEADY_PRECONDITIONED_PROBLEM_TRAITS_HPP_
#define ROM_LSPG_IMPL_STEADY_TRAITS_ROM_LSPG_STEADY_PRECONDITIONED_PROBLEM_TRAITS_HPP_

namespace pressio{ namespace rom{

//fwd declare problem class
namespace lspg{ namespace impl{ namespace steady{
template <typename ...>
class PreconditionedProblemSteady;
}}}// end namespace pressio::rom::lspg::impl::steady

namespace details{

template <
  typename fom_system_type,
  typename lspg_state_type,
  typename decoder_type,
  typename preconditioner_type
 >
struct traits<
  ::pressio::rom::lspg::impl::steady::PreconditionedProblemSteady<
    fom_system_type, lspg_state_type, decoder_type, preconditioner_type
    >
  >
{
  static const bool is_steady_lspg = true;
  static const bool is_unsteady_lspg = false;

  using common_traits =
    ::pressio::rom::lspg::impl::steady::CommonTraits<
    fom_system_type, lspg_state_type, decoder_type>;

  using fom_system_t		= typename common_traits::fom_system_t;
  using scalar_t		= typename common_traits::scalar_t;
  using fom_native_state_t	= typename common_traits::fom_native_state_t;
  using fom_state_t		= typename common_traits::fom_state_t;
  using fom_residual_t		= typename common_traits::fom_residual_t;
  using lspg_state_t		= typename common_traits::lspg_state_t;
  using lspg_native_state_t	= typename common_traits::lspg_native_state_t;
  using lspg_residual_t		= typename common_traits::lspg_residual_t;
  using decoder_t		= typename common_traits::decoder_t;
  using decoder_jac_t		= typename common_traits::decoder_jac_t;
  using lspg_jacobian_t		= typename common_traits::lspg_jacobian_t;
  using fom_state_reconstr_t	= typename common_traits::fom_state_reconstr_t;
  using fom_states_manager_t	= typename common_traits::fom_states_manager_t;
  using ud_ops_t		= typename common_traits::ud_ops_t;
  static constexpr auto binding_sentinel = common_traits::binding_sentinel;

  using preconditioner_t = preconditioner_type;
  static_assert
  (
    ::pressio::rom::lspg::constraints::steady_preconditioner<
      preconditioner_t,
      fom_state_t,
      lspg_residual_t,
      lspg_jacobian_t
      >::value,
      "Invalid preconditioner type passed to steady LSPG"
  );

  using residual_policy_t =
    ::pressio::rom::lspg::decorator::Preconditioned<
    preconditioner_t,
    ::pressio::rom::lspg::impl::steady::ResidualPolicy<
      lspg_residual_t, fom_states_manager_t>
    >;

  using jacobian_policy_t	=
    ::pressio::rom::lspg::decorator::Preconditioned<
    preconditioner_t,
    ::pressio::rom::lspg::impl::steady::JacobianPolicy<
      fom_states_manager_t, lspg_jacobian_t, decoder_t>
    >;

  using system_t = ::pressio::rom::lspg::impl::steady::System<
    scalar_t, fom_system_t, lspg_state_type,
    lspg_residual_t, lspg_jacobian_t,
    residual_policy_t, jacobian_policy_t>;
};

}}}
#endif  // ROM_LSPG_IMPL_STEADY_TRAITS_ROM_LSPG_STEADY_PRECONDITIONED_PROBLEM_TRAITS_HPP_
