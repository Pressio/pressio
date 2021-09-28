/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_traits.hpp
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

#ifndef ROM_IMPL_ROM_LSPG_UNSTEADY_TRAITS_HPP_
#define ROM_IMPL_ROM_LSPG_UNSTEADY_TRAITS_HPP_

namespace pressio{

namespace rom{ namespace lspg{ namespace impl{

//fwd declare problem class
template <int, class ...> class UnsteadyProblem;

template <class FomSystemType, class LspgStateType, class DecoderType>
struct CommonTraitsUnsteadyContTime
{
  using fom_system_type	 = FomSystemType;
  using scalar_type	 = typename fom_system_type::scalar_type;
  using lspg_state_type = LspgStateType;
  static_assert(::pressio::rom::admissible_lspg_state<lspg_state_type>::value,
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
  using fom_velocity_type = typename fom_system_type::velocity_type;
  static_assert(std::is_same<fom_state_type, fom_velocity_type>::value,
		"Currently, the fom state and residual must be of the same type");

  // ---------------------
  /* for lspg, the residual is of the same type as RHS*/
  using lspg_residual_type = fom_velocity_type;

  /* lspg_jacobian_t is type to represent dR/dx_rom where R is the residual
   * dR/dx_rom = df/dxFom * dxFom/dxRom
   *  so df/dxFom = J  and dxFom/dxRom = decoder_jacobian
   * For now, set lspg_jacobian_t to be of same type as decoder_jac_t
   * not a bad assumption since all matrices are left-applied to decoder_jac_t
   */
  using lspg_jacobian_type = decoder_jac_type;

  // fom state reconstructor
  using fom_state_reconstr_type = ::pressio::rom::FomStateReconstructor<decoder_type>;

  using fom_states_manager_type = ::pressio::rom::ManagerStencilFomStatesDynamic<
    fom_state_type, fom_state_reconstr_type, ::pressio::ode::nPlusOne>;
};


template <std::size_t num_states, class FomSystemType, class LspgStateType, class DecoderType>
struct CommonTraitsUnsteadyDiscTime
{
  using fom_system_type	 = FomSystemType;
  using scalar_type	 = typename fom_system_type::scalar_type;
  using lspg_state_type = LspgStateType;
  static_assert(::pressio::rom::admissible_lspg_state<lspg_state_type>::value,
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

  using lspg_residual_type = typename fom_system_type::discrete_time_residual_type;

  /* lspg_jacobian_t is type to represent dR/dx_rom where R is the residual
   * dR/dx_rom = df/dxFom * dxFom/dxRom
   *  so df/dxFom = J  and dxFom/dxRom = decoder_jacobian
   * For now, set lspg_jacobian_t to be of same type as decoder_jac_t
   * not a bad assumption since all matrices are left-applied to decoder_jac_t
   */
  using lspg_jacobian_type = decoder_jac_type;

  // fom state reconstructor
  using fom_state_reconstr_type = ::pressio::rom::FomStateReconstructor<decoder_type>;

  using fom_states_manager_type = ::pressio::rom::ManagerStencilFomStatesStatic<
    fom_state_type, fom_state_reconstr_type, num_states>;


};

}}} // end namespace

//=======================
// DEFAULT
//=======================

// cont-time
template <class FomSystemType, class LspgStateType, class DecoderType>
struct Traits<
  ::pressio::rom::lspg::impl::UnsteadyProblem<
    0, FomSystemType, LspgStateType, DecoderType
    >
  >
{

  using common_types = ::pressio::rom::lspg::impl::CommonTraitsUnsteadyContTime<
    FomSystemType, LspgStateType, DecoderType>;

  static constexpr auto is_cont_time = true;

  using scalar_type       = typename common_types::scalar_type;
  using fom_system_type   = typename common_types::fom_system_type;
  using fom_state_type    = typename common_types::fom_state_type;
  using fom_velocity_type = typename common_types::fom_velocity_type;

  using decoder_type		= typename common_types::decoder_type;
  using decoder_jac_type	= typename common_types::decoder_jac_type;
  using fom_state_reconstr_type = typename common_types::fom_state_reconstr_type;
  using fom_states_manager_type = typename common_types::fom_states_manager_type;

  using lspg_state_type	   = typename common_types::lspg_state_type;
  using lspg_residual_type = typename common_types::lspg_residual_type;
  using lspg_jacobian_type = typename common_types::lspg_jacobian_type;
  using size_type = typename ::pressio::Traits<lspg_state_type>::size_type;

  using residual_policy_type =
    ::pressio::rom::lspg::impl::UnsteadyResidualPolicy<
    lspg_residual_type, fom_states_manager_type, fom_system_type
    >;

  using jacobian_policy_type =
    ::pressio::rom::lspg::impl::UnsteadyJacobianPolicy<
    lspg_jacobian_type, fom_states_manager_type, decoder_type, fom_system_type
    >;

  using stepper_type = typename ::pressio::ode::impl::ImplicitCompose<
    lspg_state_type, residual_policy_type &, jacobian_policy_type &>::type;
};

// discrete-time
template <
  std::size_t num_states,
  class FomSystemType,
  class LspgStateType,
  class DecoderType
  >
struct Traits<
  ::pressio::rom::lspg::impl::UnsteadyProblem<
    1, FomSystemType, LspgStateType, DecoderType,
    ::pressio::ode::StepperTotalNumberOfStates<num_states>
    >
  >
{

  using common_types = ::pressio::rom::lspg::impl::CommonTraitsUnsteadyDiscTime<
    num_states, FomSystemType, LspgStateType, DecoderType>;

  static constexpr auto is_cont_time = false;

  using scalar_type       = typename common_types::scalar_type;
  using fom_system_type   = typename common_types::fom_system_type;
  using fom_state_type    = typename common_types::fom_state_type;

  using decoder_type		= typename common_types::decoder_type;
  using decoder_jac_type	= typename common_types::decoder_jac_type;
  using fom_state_reconstr_type = typename common_types::fom_state_reconstr_type;
  using fom_states_manager_type = typename common_types::fom_states_manager_type;

  using lspg_state_type	   = typename common_types::lspg_state_type;
  using lspg_residual_type = typename common_types::lspg_residual_type;
  using lspg_jacobian_type = typename common_types::lspg_jacobian_type;
  using size_type = typename ::pressio::Traits<lspg_state_type>::size_type;

  using rom_system_type = ::pressio::rom::lspg::impl::DiscreteTimeDefaultSystem<
    scalar_type, lspg_state_type, lspg_residual_type, lspg_jacobian_type,
    fom_system_type, decoder_type, fom_states_manager_type>;

  using stepper_type = typename ::pressio::ode::impl::ImplicitComposeArb<
    num_states, const rom_system_type &, lspg_state_type>::type;
};


//=======================
// PRECONDITIONED DEFAULT
//=======================

// cont-time
template <
  class FomSystemType,
  class LspgStateType,
  class DecoderType,
  class PreconditionerType
  >
struct Traits<
  ::pressio::rom::lspg::impl::UnsteadyProblem<
    2, FomSystemType, LspgStateType, DecoderType, PreconditionerType
    >
  >
{

  using common_types = ::pressio::rom::lspg::impl::CommonTraitsUnsteadyContTime<
    FomSystemType, LspgStateType, DecoderType>;

  static constexpr auto is_cont_time = true;

  using scalar_type       = typename common_types::scalar_type;
  using fom_system_type   = typename common_types::fom_system_type;
  using fom_state_type    = typename common_types::fom_state_type;
  using fom_velocity_type = typename common_types::fom_velocity_type;

  using decoder_type		= typename common_types::decoder_type;
  using decoder_jac_type	= typename common_types::decoder_jac_type;
  using fom_state_reconstr_type = typename common_types::fom_state_reconstr_type;
  using fom_states_manager_type = typename common_types::fom_states_manager_type;

  using lspg_state_type	   = typename common_types::lspg_state_type;
  using lspg_residual_type = typename common_types::lspg_residual_type;
  using lspg_jacobian_type = typename common_types::lspg_jacobian_type;
  using size_type = typename ::pressio::Traits<lspg_state_type>::size_type;

  using preconditioner_type = PreconditionerType;

  using residual_policy_type =
    ::pressio::rom::lspg::impl::PrecDecoratorResidual<
    preconditioner_type,
    ::pressio::rom::lspg::impl::UnsteadyResidualPolicy<
      lspg_residual_type, fom_states_manager_type, fom_system_type
      >
    >;

  using jacobian_policy_type =
    ::pressio::rom::lspg::impl::PrecDecoratorJacobian<
    preconditioner_type,
    ::pressio::rom::lspg::impl::UnsteadyJacobianPolicy<
      lspg_jacobian_type, fom_states_manager_type, decoder_type, fom_system_type
      >
    >;

  using stepper_type = typename ::pressio::ode::impl::ImplicitCompose<
    lspg_state_type, residual_policy_type &, jacobian_policy_type &>::type;
};

// discrete-time
template <
  std::size_t num_states,
  class FomSystemType,
  class LspgStateType,
  class DecoderType,
  class PreconditionerType
  >
struct Traits<
  ::pressio::rom::lspg::impl::UnsteadyProblem<
    3, FomSystemType, LspgStateType, DecoderType,
    ::pressio::ode::StepperTotalNumberOfStates<num_states>, PreconditionerType
    >
  >
{

  using common_types = ::pressio::rom::lspg::impl::CommonTraitsUnsteadyDiscTime<
    num_states, FomSystemType, LspgStateType, DecoderType>;

  static constexpr auto is_cont_time = false;

  using scalar_type       = typename common_types::scalar_type;
  using fom_system_type   = typename common_types::fom_system_type;
  using fom_state_type    = typename common_types::fom_state_type;

  using decoder_type		= typename common_types::decoder_type;
  using decoder_jac_type	= typename common_types::decoder_jac_type;
  using fom_state_reconstr_type = typename common_types::fom_state_reconstr_type;
  using fom_states_manager_type = typename common_types::fom_states_manager_type;

  using lspg_state_type	   = typename common_types::lspg_state_type;
  using lspg_residual_type = typename common_types::lspg_residual_type;
  using lspg_jacobian_type = typename common_types::lspg_jacobian_type;
  using size_type = typename ::pressio::Traits<lspg_state_type>::size_type;

  using preconditioner_type = PreconditionerType;

  using rom_system_type =
    ::pressio::rom::lspg::impl::PrecDecoratorDiscreteTimeSystem<
    preconditioner_type,
    ::pressio::rom::lspg::impl::DiscreteTimeDefaultSystem<
      scalar_type, lspg_state_type, lspg_residual_type, lspg_jacobian_type,
      fom_system_type, decoder_type, fom_states_manager_type
    >
    >;

  using stepper_type = typename ::pressio::ode::impl::ImplicitComposeArb<
    num_states, const rom_system_type &, lspg_state_type>::type;
};


//=======================
// MASKED
//=======================

// cont-time
template <
  class FomSystemType,
  class LspgStateType,
  class DecoderType,
  class MaskerType
  >
struct Traits<
  ::pressio::rom::lspg::impl::UnsteadyProblem<
    4, FomSystemType, LspgStateType, DecoderType, MaskerType
    >
  >
{

  using common_types = ::pressio::rom::lspg::impl::CommonTraitsUnsteadyContTime<
    FomSystemType, LspgStateType, DecoderType>;

  static constexpr auto is_cont_time = true;

  using scalar_type       = typename common_types::scalar_type;
  using fom_system_type   = typename common_types::fom_system_type;
  using fom_state_type    = typename common_types::fom_state_type;
  using fom_velocity_type = typename common_types::fom_velocity_type;

  using decoder_type		= typename common_types::decoder_type;
  using decoder_jac_type	= typename common_types::decoder_jac_type;
  using fom_state_reconstr_type = typename common_types::fom_state_reconstr_type;
  using fom_states_manager_type = typename common_types::fom_states_manager_type;

  using lspg_state_type	   = typename common_types::lspg_state_type;
  using lspg_residual_type = typename common_types::lspg_residual_type;
  using lspg_jacobian_type = typename common_types::lspg_jacobian_type;
  using size_type = typename ::pressio::Traits<lspg_state_type>::size_type;

  using masker_type = MaskerType;

  using residual_policy_type =
    ::pressio::rom::lspg::impl::MaskDecoratorResidual<
    masker_type,
    ::pressio::rom::lspg::impl::UnsteadyResidualPolicy<
      lspg_residual_type, fom_states_manager_type, fom_system_type
      >
    >;

  using jacobian_policy_type =
    ::pressio::rom::lspg::impl::MaskDecoratorJacobian<
    masker_type,
    ::pressio::rom::lspg::impl::UnsteadyJacobianPolicy<
      lspg_jacobian_type, fom_states_manager_type, decoder_type, fom_system_type
      >
    >;

  using stepper_type = typename ::pressio::ode::impl::ImplicitCompose<
    lspg_state_type, residual_policy_type &, jacobian_policy_type &>::type;
};

// discrete-time
template <
  std::size_t num_states,
  class FomSystemType,
  class LspgStateType,
  class DecoderType,
  class MaskerType
  >
struct Traits<
  ::pressio::rom::lspg::impl::UnsteadyProblem<
    5, FomSystemType, LspgStateType, DecoderType,
    ::pressio::ode::StepperTotalNumberOfStates<num_states>, MaskerType
    >
  >
{

  using common_types = ::pressio::rom::lspg::impl::CommonTraitsUnsteadyDiscTime<
    num_states, FomSystemType, LspgStateType, DecoderType>;

  static constexpr auto is_cont_time = false;

  using scalar_type       = typename common_types::scalar_type;
  using fom_system_type   = typename common_types::fom_system_type;
  using fom_state_type    = typename common_types::fom_state_type;

  using decoder_type		= typename common_types::decoder_type;
  using decoder_jac_type	= typename common_types::decoder_jac_type;
  using fom_state_reconstr_type = typename common_types::fom_state_reconstr_type;
  using fom_states_manager_type = typename common_types::fom_states_manager_type;

  using lspg_state_type	   = typename common_types::lspg_state_type;
  using lspg_residual_type = typename common_types::lspg_residual_type;
  using lspg_jacobian_type = typename common_types::lspg_jacobian_type;
  using size_type = typename ::pressio::Traits<lspg_state_type>::size_type;

  using masker_type = MaskerType;

  using rom_system_type =
    ::pressio::rom::lspg::impl::MaskDecoratorDiscreteTimeSystem<
    masker_type,
    ::pressio::rom::lspg::impl::DiscreteTimeDefaultSystem<
      scalar_type, lspg_state_type, lspg_residual_type, lspg_jacobian_type,
      fom_system_type, decoder_type, fom_states_manager_type
    >
    >;

  using stepper_type = typename ::pressio::ode::impl::ImplicitComposeArb<
    num_states, const rom_system_type &, lspg_state_type>::type;
};

//=======================
// PRECONDITIONED MASKED
//=======================

// cont-time
template <
  class FomSystemType,
  class LspgStateType,
  class DecoderType,
  class MaskerType,
  class PreconditionerType
  >
struct Traits<
  ::pressio::rom::lspg::impl::UnsteadyProblem<
    6, FomSystemType, LspgStateType, DecoderType, MaskerType, PreconditionerType
    >
  >
{

  using common_types = ::pressio::rom::lspg::impl::CommonTraitsUnsteadyContTime<
    FomSystemType, LspgStateType, DecoderType>;

  static constexpr auto is_cont_time = true;

  using scalar_type       = typename common_types::scalar_type;
  using fom_system_type   = typename common_types::fom_system_type;
  using fom_state_type    = typename common_types::fom_state_type;
  using fom_velocity_type = typename common_types::fom_velocity_type;

  using decoder_type		= typename common_types::decoder_type;
  using decoder_jac_type	= typename common_types::decoder_jac_type;
  using fom_state_reconstr_type = typename common_types::fom_state_reconstr_type;
  using fom_states_manager_type = typename common_types::fom_states_manager_type;

  using lspg_state_type	   = typename common_types::lspg_state_type;
  using lspg_residual_type = typename common_types::lspg_residual_type;
  using lspg_jacobian_type = typename common_types::lspg_jacobian_type;
  using size_type = typename ::pressio::Traits<lspg_state_type>::size_type;

  using preconditioner_type = PreconditionerType;
  using masker_type = MaskerType;

  using residual_policy_type =
    ::pressio::rom::lspg::impl::PrecDecoratorResidual<
    preconditioner_type,
    ::pressio::rom::lspg::impl::MaskDecoratorResidual<
      masker_type,
      ::pressio::rom::lspg::impl::UnsteadyResidualPolicy<
	lspg_residual_type, fom_states_manager_type, fom_system_type
	>
      >
    >;

  using jacobian_policy_type =
    ::pressio::rom::lspg::impl::PrecDecoratorJacobian<
    preconditioner_type,
    ::pressio::rom::lspg::impl::MaskDecoratorJacobian<
      masker_type,
      ::pressio::rom::lspg::impl::UnsteadyJacobianPolicy<
	lspg_jacobian_type, fom_states_manager_type, decoder_type, fom_system_type
	>
      >
    >;

  using stepper_type = typename ::pressio::ode::impl::ImplicitCompose<
    lspg_state_type, residual_policy_type &, jacobian_policy_type &>::type;
};

// discrete-time
template <
  std::size_t num_states,
  class FomSystemType,
  class LspgStateType,
  class DecoderType,
  class MaskerType,
  class PreconditionerType
  >
struct Traits<
  ::pressio::rom::lspg::impl::UnsteadyProblem<
    7, FomSystemType, LspgStateType, DecoderType,
    ::pressio::ode::StepperTotalNumberOfStates<num_states>,
    MaskerType, PreconditionerType
    >
  >
{

  using common_types = ::pressio::rom::lspg::impl::CommonTraitsUnsteadyDiscTime<
    num_states, FomSystemType, LspgStateType, DecoderType>;

  static constexpr auto is_cont_time = false;

  using scalar_type       = typename common_types::scalar_type;
  using fom_system_type   = typename common_types::fom_system_type;
  using fom_state_type    = typename common_types::fom_state_type;

  using decoder_type		= typename common_types::decoder_type;
  using decoder_jac_type	= typename common_types::decoder_jac_type;
  using fom_state_reconstr_type = typename common_types::fom_state_reconstr_type;
  using fom_states_manager_type = typename common_types::fom_states_manager_type;

  using lspg_state_type	   = typename common_types::lspg_state_type;
  using lspg_residual_type = typename common_types::lspg_residual_type;
  using lspg_jacobian_type = typename common_types::lspg_jacobian_type;
  using size_type = typename ::pressio::Traits<lspg_state_type>::size_type;

  using preconditioner_type = PreconditionerType;
  using masker_type = MaskerType;

  using rom_system_type =
    ::pressio::rom::lspg::impl::PrecDecoratorDiscreteTimeSystem<
    preconditioner_type,
    ::pressio::rom::lspg::impl::MaskDecoratorDiscreteTimeSystem<
      masker_type,
      ::pressio::rom::lspg::impl::DiscreteTimeDefaultSystem<
	scalar_type, lspg_state_type, lspg_residual_type, lspg_jacobian_type,
	fom_system_type, decoder_type, fom_states_manager_type
	>
      >
    >;

  using stepper_type = typename ::pressio::ode::impl::ImplicitComposeArb<
    num_states, const rom_system_type &, lspg_state_type>::type;
};


//=======================
// HYPER-REDUCED
//=======================

// cont-time
template <
  class FomSystemType,
  class LspgStateType,
  class DecoderType,
  class HypRedOperatorUpdater
  >
struct Traits<
  ::pressio::rom::lspg::impl::UnsteadyProblem<
    8, FomSystemType, LspgStateType, DecoderType, HypRedOperatorUpdater
    >
  >
{

  using common_types = ::pressio::rom::lspg::impl::CommonTraitsUnsteadyContTime<
    FomSystemType, LspgStateType, DecoderType>;

  static constexpr auto is_cont_time = true;

  using scalar_type       = typename common_types::scalar_type;
  using fom_system_type   = typename common_types::fom_system_type;
  using fom_state_type    = typename common_types::fom_state_type;
  using fom_velocity_type = typename common_types::fom_velocity_type;

  using decoder_type		= typename common_types::decoder_type;
  using decoder_jac_type	= typename common_types::decoder_jac_type;
  using fom_state_reconstr_type = typename common_types::fom_state_reconstr_type;
  using fom_states_manager_type = typename common_types::fom_states_manager_type;

  using lspg_state_type	   = typename common_types::lspg_state_type;
  using lspg_residual_type = typename common_types::lspg_residual_type;
  using lspg_jacobian_type = typename common_types::lspg_jacobian_type;
  using size_type = typename ::pressio::Traits<lspg_state_type>::size_type;

  using hyperreduced_operator_type = HypRedOperatorUpdater;

  using residual_policy_type =
    ::pressio::rom::lspg::impl::UnsteadyHypRedResidualPolicy<
    lspg_residual_type, fom_states_manager_type,
    fom_system_type, hyperreduced_operator_type
    >;

  using jacobian_policy_type =
    ::pressio::rom::lspg::impl::UnsteadyHypRedJacobianPolicy<
    lspg_jacobian_type, fom_states_manager_type, decoder_type,
    fom_system_type, hyperreduced_operator_type
    >;

  using stepper_type = typename ::pressio::ode::impl::ImplicitCompose<
    lspg_state_type, residual_policy_type &, jacobian_policy_type &>::type;
};


//===============================
// PRECONDITIONED HYPER-REDUCED
//===============================

// cont-time
template <
  class FomSystemType,
  class LspgStateType,
  class DecoderType,
  class HypRedOperatorUpdater,
  class PreconditionerType
  >
struct Traits<
  ::pressio::rom::lspg::impl::UnsteadyProblem<
    9, FomSystemType, LspgStateType, DecoderType, HypRedOperatorUpdater, PreconditionerType
    >
  >
{

  using common_types = ::pressio::rom::lspg::impl::CommonTraitsUnsteadyContTime<
    FomSystemType, LspgStateType, DecoderType>;

  static constexpr auto binding_sentinel = common_types::binding_sentinel;
  static constexpr auto is_cont_time = true;

  using scalar_type       = typename common_types::scalar_type;
  using fom_system_type   = typename common_types::fom_system_type;
  using fom_state_type    = typename common_types::fom_state_type;
  using fom_velocity_type = typename common_types::fom_velocity_type;

  using decoder_type		= typename common_types::decoder_type;
  using decoder_jac_type	= typename common_types::decoder_jac_type;
  using fom_state_reconstr_type = typename common_types::fom_state_reconstr_type;
  using fom_states_manager_type = typename common_types::fom_states_manager_type;

  using lspg_state_type	   = typename common_types::lspg_state_type;
  using lspg_residual_type = typename common_types::lspg_residual_type;
  using lspg_jacobian_type = typename common_types::lspg_jacobian_type;
  using size_type = typename ::pressio::Traits<lspg_state_type>::size_type;

  using hyperreduced_operator_type = HypRedOperatorUpdater;
  using preconditioner_type = PreconditionerType;

  using residual_policy_type =
    ::pressio::rom::lspg::impl::PrecDecoratorResidual<
    preconditioner_type,
    ::pressio::rom::lspg::impl::UnsteadyHypRedResidualPolicy<
      lspg_residual_type, fom_states_manager_type,
      fom_system_type, hyperreduced_operator_type
      >
    >;

  using jacobian_policy_type =
    ::pressio::rom::lspg::impl::PrecDecoratorJacobian<
    preconditioner_type,
    ::pressio::rom::lspg::impl::UnsteadyHypRedJacobianPolicy<
      lspg_jacobian_type, fom_states_manager_type, decoder_type,
      fom_system_type, hyperreduced_operator_type
      >
    >;

  using stepper_type = typename ::pressio::ode::impl::ImplicitCompose<
    lspg_state_type, residual_policy_type &, jacobian_policy_type &>::type;
};

}//end  namespace pressio
#endif  // ROM_IMPL_ROM_LSPG_UNSTEADY_TRAITS_HPP_
