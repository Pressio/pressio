/*
//@HEADER
// ************************************************************************
//
// rom_lspg_steady_traits.hpp
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

#ifndef ROM_IMPL_ROM_LSPG_STEADY_TRAITS_HPP_
#define ROM_IMPL_ROM_LSPG_STEADY_TRAITS_HPP_

namespace pressio{

namespace rom{ namespace lspg{ namespace impl{

//fwd declare problem class
template <int, class ...> class SteadyProblem;

template <class FomSystemType, class LspgStateType, class DecoderType>
struct CommonTraitsSteady
{
  using fom_system_type	 = FomSystemType;
  using scalar_type	 = typename fom_system_type::scalar_type;
  using lspg_state_type = LspgStateType;

  static_assert
  (::pressio::rom::admissible_lspg_state<lspg_state_type>::value,
   "Invalid lspg state type");

  using decoder_type = DecoderType;
  using decoder_jac_type = typename decoder_type::jacobian_type;

  // ---------------------
  // detect fom state type from decoder
  // ensure it is consistent with the fom_state_type from the app
  using fom_state_from_decoder_type = typename decoder_type::fom_state_type;
  using fom_state_from_adapter_type = typename fom_system_type::state_type;
  static_assert
  (std::is_same<fom_state_from_decoder_type, fom_state_from_adapter_type>::value,
   "The fom state type detected from the fom adapter must match the fom state type used in the decoder");
  using fom_state_type = fom_state_from_decoder_type;

  // ---------------------
  using fom_residual_type = fom_state_type;
  using fom_residual_from_adapter_type = typename fom_system_type::residual_type;
  static_assert
  (std::is_same<fom_state_type, fom_residual_from_adapter_type>::value,
   "Currently, the fom state and residual must be of the same type");

  // ---------------------
  /* for steady lspg, the lspg residual is literally the rhs (aka velocity fo us) */
  using lspg_residual_type = fom_residual_type;

  /* lspg_jacobian_t is type to represent dR/dx_rom where R is the residual
   * dR/dx_rom = df/dxFom * dxFom/dxRom
   *  so df/dxFom = J  and dxFom/dxRom = decoder_jacobian
   * For now, set lspg_jacobian_t to be of same type as decoder_jac_t
   * not a bad assumption since all matrices are left-applied to decoder_jac_t
   */
  using lspg_jacobian_type = decoder_jac_type;

  // ---------------------
  // fom state reconstructor
  using fom_state_reconstr_type = ::pressio::rom::FomStateReconstructor<decoder_type>;

  // for steady lspg, we only need to store one FOM state
  using fom_states_manager_type = ::pressio::rom::ManagerFomStatesSteady<
    fom_state_type, fom_state_reconstr_type>;

  static constexpr bool binding_sentinel = false;
};

}}} // end namespace pressio::rom::galekin::impl

//=======================
// DEFAULT
//=======================
template <
  class FomSystemType,
  class LspgStateType,
  class DecoderType
  >
struct Traits<
  ::pressio::rom::lspg::impl::SteadyProblem<
    0, FomSystemType, LspgStateType, DecoderType
    >
  >
{
  using common_types = ::pressio::rom::lspg::impl::CommonTraitsSteady<
    FomSystemType, LspgStateType, DecoderType>;

  static constexpr auto binding_sentinel = common_types::binding_sentinel;
  using fom_system_type    = typename common_types::fom_system_type;
  using scalar_type    = typename common_types::scalar_type;

  using fom_state_type   = typename common_types::fom_state_type;
  using fom_residual_type    = typename common_types::fom_residual_type;

  using decoder_type   = typename common_types::decoder_type;
  using decoder_jac_type   = typename common_types::decoder_jac_type;
  using fom_state_reconstr_type  = typename common_types::fom_state_reconstr_type;
  using fom_states_manager_type  = typename common_types::fom_states_manager_type;

  using lspg_state_type  = typename common_types::lspg_state_type;
  using lspg_residual_type = typename common_types::lspg_residual_type;
  using lspg_jacobian_type = typename common_types::lspg_jacobian_type;
  using size_type = typename ::pressio::Traits<lspg_state_type>::size_type;

  using residual_policy_type =
    ::pressio::rom::lspg::impl::SteadyResidualPolicy<
    fom_states_manager_type, lspg_residual_type, fom_system_type
    >;

  using jacobian_policy_type =
    ::pressio::rom::lspg::impl::SteadyJacobianPolicy<
    fom_states_manager_type, decoder_jac_type, decoder_type,fom_system_type
    >;

  using steady_system_type = ::pressio::rom::lspg::impl::SteadySystem<
    scalar_type, lspg_state_type, lspg_residual_type, lspg_jacobian_type,
    residual_policy_type, jacobian_policy_type>;
};

//=======================
// MASKED
//=======================
template <
  class FomSystemType,
  class LspgStateType,
  class DecoderType,
  class MaskerType
  >
struct Traits<
  ::pressio::rom::lspg::impl::SteadyProblem<
    1, FomSystemType, LspgStateType, DecoderType, MaskerType
    >
  >
{
  using common_types = ::pressio::rom::lspg::impl::CommonTraitsSteady<
    FomSystemType, LspgStateType, DecoderType>;

  static constexpr auto binding_sentinel = common_types::binding_sentinel;
  using fom_system_type    = typename common_types::fom_system_type;
  using scalar_type    = typename common_types::scalar_type;

  using fom_state_type   = typename common_types::fom_state_type;
  using fom_residual_type    = typename common_types::fom_residual_type;

  using decoder_type   = typename common_types::decoder_type;
  using decoder_jac_type   = typename common_types::decoder_jac_type;
  using fom_state_reconstr_type  = typename common_types::fom_state_reconstr_type;
  using fom_states_manager_type  = typename common_types::fom_states_manager_type;

  using lspg_state_type  = typename common_types::lspg_state_type;
  using lspg_residual_type = typename common_types::lspg_residual_type;
  using lspg_jacobian_type = typename common_types::lspg_jacobian_type;
  using size_type = typename ::pressio::Traits<lspg_state_type>::size_type;

  using masker_type = MaskerType;
  static_assert
  (::pressio::rom::lspg_steady_masker<
   masker_type, fom_residual_type, lspg_jacobian_type>::value,
   "Invalid masker passed for steady LSPG");

  using residual_policy_type =
    ::pressio::rom::lspg::impl::MaskDecoratorSteady<
    lspg_residual_type, masker_type,
    ::pressio::rom::lspg::impl::SteadyResidualPolicy<
      fom_states_manager_type, lspg_residual_type, fom_system_type
      >
    >;

  using jacobian_policy_type =
    ::pressio::rom::lspg::impl::MaskDecoratorSteady<
    lspg_jacobian_type, masker_type,
    ::pressio::rom::lspg::impl::SteadyJacobianPolicy<
      fom_states_manager_type, lspg_jacobian_type, decoder_type,fom_system_type
      >
    >;

  using steady_system_type = ::pressio::rom::lspg::impl::SteadySystem<
    scalar_type, lspg_state_type, lspg_residual_type, lspg_jacobian_type,
    residual_policy_type, jacobian_policy_type>;
};

//=======================
// PRECONDITIONED DEFAULT
//=======================
template <
  class FomSystemType,
  class LspgStateType,
  class DecoderType,
  class PreconditionerType
  >
struct Traits<
  ::pressio::rom::lspg::impl::SteadyProblem<
    2, FomSystemType, LspgStateType, DecoderType, PreconditionerType
    >
  >
{
  using common_types = ::pressio::rom::lspg::impl::CommonTraitsSteady<
    FomSystemType, LspgStateType, DecoderType>;

  static constexpr auto binding_sentinel = common_types::binding_sentinel;
  using fom_system_type    = typename common_types::fom_system_type;
  using scalar_type    = typename common_types::scalar_type;

  using fom_state_type   = typename common_types::fom_state_type;
  using fom_residual_type    = typename common_types::fom_residual_type;

  using decoder_type   = typename common_types::decoder_type;
  using decoder_jac_type   = typename common_types::decoder_jac_type;
  using fom_state_reconstr_type  = typename common_types::fom_state_reconstr_type;
  using fom_states_manager_type  = typename common_types::fom_states_manager_type;

  using lspg_state_type  = typename common_types::lspg_state_type;
  using lspg_residual_type = typename common_types::lspg_residual_type;
  using lspg_jacobian_type = typename common_types::lspg_jacobian_type;
  using size_type = typename ::pressio::Traits<lspg_state_type>::size_type;

  using preconditioner_type = PreconditionerType;
  static_assert
  (::pressio::rom::lspg_steady_preconditioner<
   preconditioner_type, fom_state_type, lspg_residual_type, lspg_jacobian_type>::value,
   "Invalid masker passed for steady LSPG");

  using residual_policy_type =
    ::pressio::rom::lspg::impl::PreconditionDecoratorSteady<
    lspg_residual_type, preconditioner_type,
    ::pressio::rom::lspg::impl::SteadyResidualPolicy<
      fom_states_manager_type, lspg_residual_type, fom_system_type
      >
    >;

  using jacobian_policy_type = ::pressio::rom::lspg::impl::PreconditionDecoratorSteady<
    lspg_jacobian_type, preconditioner_type,
    ::pressio::rom::lspg::impl::SteadyJacobianPolicy<
      fom_states_manager_type, lspg_jacobian_type, decoder_type,fom_system_type
      >
    >;

  using steady_system_type = ::pressio::rom::lspg::impl::SteadySystem<
    scalar_type, lspg_state_type, lspg_residual_type, lspg_jacobian_type,
    residual_policy_type, jacobian_policy_type>;
};

//=======================
// PRECONDITIONED MASKED
//=======================
template <
  class FomSystemType,
  class LspgStateType,
  class DecoderType,
  class MaskerType,
  class PreconditionerType
  >
struct Traits<
  ::pressio::rom::lspg::impl::SteadyProblem<
    3, FomSystemType, LspgStateType, DecoderType, MaskerType, PreconditionerType
    >
  >
{
  using common_types = ::pressio::rom::lspg::impl::CommonTraitsSteady<
    FomSystemType, LspgStateType, DecoderType>;

  static constexpr auto binding_sentinel = common_types::binding_sentinel;
  using fom_system_type    = typename common_types::fom_system_type;
  using scalar_type    = typename common_types::scalar_type;

  using fom_state_type   = typename common_types::fom_state_type;
  using fom_residual_type    = typename common_types::fom_residual_type;

  using decoder_type   = typename common_types::decoder_type;
  using decoder_jac_type   = typename common_types::decoder_jac_type;
  using fom_state_reconstr_type  = typename common_types::fom_state_reconstr_type;
  using fom_states_manager_type  = typename common_types::fom_states_manager_type;

  using lspg_state_type  = typename common_types::lspg_state_type;
  using lspg_residual_type = typename common_types::lspg_residual_type;
  using lspg_jacobian_type = typename common_types::lspg_jacobian_type;
  using size_type = typename ::pressio::Traits<lspg_state_type>::size_type;

  using masker_type = MaskerType;
  static_assert
  (::pressio::rom::lspg_steady_masker<
   masker_type, fom_residual_type, lspg_jacobian_type>::value,
   "Invalid masker passed for steady LSPG");

  using preconditioner_type = PreconditionerType;
  static_assert
  (::pressio::rom::lspg_steady_preconditioner<
   preconditioner_type, fom_state_type, lspg_residual_type, lspg_jacobian_type>::value,
   "Invalid masker passed for steady LSPG");

  using residual_policy_type =
    ::pressio::rom::lspg::impl::PreconditionDecoratorSteady<
    lspg_residual_type, preconditioner_type,
    ::pressio::rom::lspg::impl::MaskDecoratorSteady<
      lspg_residual_type, masker_type,
      ::pressio::rom::lspg::impl::SteadyResidualPolicy<
	fom_states_manager_type, lspg_residual_type, fom_system_type
	>
      >
    >;

  using jacobian_policy_type =
    ::pressio::rom::lspg::impl::PreconditionDecoratorSteady<
    lspg_jacobian_type, preconditioner_type,
    ::pressio::rom::lspg::impl::MaskDecoratorSteady<
      lspg_jacobian_type, masker_type,
      ::pressio::rom::lspg::impl::SteadyJacobianPolicy<
	fom_states_manager_type, lspg_jacobian_type, decoder_type,fom_system_type
	>
      >
    >;

  using steady_system_type = ::pressio::rom::lspg::impl::SteadySystem<
    scalar_type, lspg_state_type, lspg_residual_type, lspg_jacobian_type,
    residual_policy_type, jacobian_policy_type>;
};

}//end  namespace pressio
#endif  // ROM_IMPL_ROM_LSPG_STEADY_TRAITS_HPP_
