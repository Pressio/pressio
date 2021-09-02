/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_common_traits.hpp
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

#ifndef ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_TRAITS_ROM_GALERKIN_COMMON_TRAITS_HPP_
#define ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_TRAITS_ROM_GALERKIN_COMMON_TRAITS_HPP_

namespace pressio{

namespace rom{ namespace galerkin{ namespace impl{
//fwd declare problem class
template <int, class ...> class Problem;

template<typename tag>
struct _num_fom_states_needed{
  static constexpr std::size_t value = 1;
};

template<>
struct _num_fom_states_needed<::pressio::ode::BDF1>{
  static constexpr std::size_t value = 1;
};

template<>
struct _num_fom_states_needed<::pressio::ode::BDF2>{
  static constexpr std::size_t value = 1;
};

template<>
struct _num_fom_states_needed<::pressio::ode::CrankNicolson>{
  static constexpr std::size_t value = 2;
};

template <
  class OdeTag,
  class FomSystemType,
  class GalerkinStateType,
  class DecoderType
  >
struct CommonTraitsContinuousTimeApi
{
  using fom_system_type	 = FomSystemType;
  using scalar_type	 = typename fom_system_type::scalar_type;
  using galerkin_state_type = GalerkinStateType;

  static_assert
  (::pressio::rom::admissible_galerkin_state<galerkin_state_type>::value,
   "Invalid galerkin state type");

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
  // we don't allow fom state and fom velocity to have different types
  // but need to make sure this assumption is consistent with fom class
  using fom_velocity_type = fom_state_type;
  using fom_velocity_from_adapter_type = typename fom_system_type::velocity_type;
  static_assert
  (std::is_same<fom_state_type, fom_velocity_from_adapter_type>::value,
   "Currently, the fom state and velocity must be of the same type");

  // ---------------------
  // fom state reconstructor
  using fom_state_reconstr_type = ::pressio::rom::FomStateReconstructor<decoder_type>;

  // ---------------------------------------------------------------
  /* fom states manager

     - This is Galerkin, so the FOM states are only needed when we query the FOM velocity.

     - For explicit time stepping, we only need to store one fom state that
     we reconstruct every time we need to compute the FOM velocity

     - For implicit time stepping, we might need to store more than one
     depending on how many FOM velocity evaluations are needed.
     BDF1 and BDF2 need FOM velocity at n+1, so one FOM state.
     CrankNicolson needs two evaluations of the FOM velocity
     at n and n+1, so we need two FOM states.
   */
  static constexpr auto nstates = _num_fom_states_needed<OdeTag>::value;
  using fom_states_manager_type = mpl::conditional_t<
    ::pressio::ode::is_explicit_stepper_tag<OdeTag>::value,
    ::pressio::rom::ManagerFomStatesUnsteadyExplicit<
      fom_state_type, fom_state_reconstr_type, nstates>,
    ::pressio::rom::ManagerFomStatesUnsteadyImplicit<
      fom_state_type, fom_state_reconstr_type, nstates>
    >;

  // ---------------------
  // sentinel to tell if we are doing bindings for p4py:
  // always false if pybind is disabled, otherwise detect from galerkin state
  static constexpr bool binding_sentinel = false;
// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
//     ::pressio::containers::predicates::is_tensor_wrapper_pybind<galerkin_state_t>::value;
// #else
//   false;
// #endif
};


template <
  std::size_t num_states,
  class FomSystemType,
  class GalerkinStateType,
  class DecoderType
  >
struct CommonTraitsDiscreteTimeApi
{
  using fom_system_type	 = FomSystemType;
  using scalar_type	 = typename fom_system_type::scalar_type;
  using galerkin_state_type = GalerkinStateType;

  static_assert
  (::pressio::rom::admissible_galerkin_state<galerkin_state_type>::value,
   "Invalid galerkin state type");

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
  using fom_discrete_time_residual_type = typename fom_system_type::discrete_time_residual_type;

  // ---------------------
  // fom state reconstructor
  using fom_state_reconstr_type = ::pressio::rom::FomStateReconstructor<decoder_type>;

  static constexpr auto nstates = num_states;
  using fom_states_manager_type = ::pressio::rom::ManagerFomStatesUnsteadyImplicit<
    fom_state_type, fom_state_reconstr_type, nstates>;

  static constexpr bool binding_sentinel = false;
};

}}} // end namespace pressio::rom::galekin::impl



//=======================
//
// DEFAULT
//
//=======================

// conttime explicit
template <
  class StepperTag,
  class FomSystemType,
  class GalerkinStateType,
  class DecoderType
  >
struct Traits<
  ::pressio::rom::galerkin::impl::Problem<
    0, StepperTag, FomSystemType, GalerkinStateType, DecoderType
    >
  >
{
  using common_types = ::pressio::rom::galerkin::impl::CommonTraitsContinuousTimeApi<
    StepperTag, FomSystemType, GalerkinStateType, DecoderType>;
  using fom_system_type		= typename common_types::fom_system_type;
  using scalar_type		= typename common_types::scalar_type;
  using fom_state_type		= typename common_types::fom_state_type;
  using fom_velocity_type		= typename common_types::fom_velocity_type;
  using galerkin_state_type	= typename common_types::galerkin_state_type;
  // for now, the velocity type is same as state
  using galerkin_velocity_type = galerkin_state_type;

  using decoder_type		= typename common_types::decoder_type;
  using decoder_jac_type		= typename common_types::decoder_jac_type;
  using fom_state_reconstr_type	= typename common_types::fom_state_reconstr_type;
  using fom_states_manager_type	= typename common_types::fom_states_manager_type;
  static constexpr auto binding_sentinel = common_types::binding_sentinel;
  using size_type = typename ::pressio::Traits<galerkin_state_type>::size_type;

  // default galerkin uses default projector, so decoderJac^T
  using projector_type = ::pressio::rom::galerkin::impl::DefaultProjector<decoder_type>;

  using rom_system_type =
    ::pressio::rom::galerkin::impl::VelocityOnlySystemUsing<
    scalar_type, galerkin_state_type, galerkin_velocity_type,
    ::pressio::rom::galerkin::impl::Projected<
      projector_type, size_type, 1,
      ::pressio::rom::galerkin::impl::DefaultFomVelocityEvaluator<
        fom_states_manager_type, fom_velocity_type, fom_system_type>
      >
    >;

  // two things to fix here:
  // 1. use the public API not impl
  // 2. we need to use const rom_sytem_type & because the ROM problem
  // owns the system and the stepper only references it
  using stepper_type = typename ::pressio::ode::impl::ExplicitCompose<
    StepperTag, galerkin_state_type, const rom_system_type &>::type;
};

// conttime implicit
template <
  class StepperTag,
  class FomSystemType,
  class GalerkinStateType,
  class GalerkinResidualType,
  class GalerkinJacobianType,
  class DecoderType
  >
struct Traits<
  ::pressio::rom::galerkin::impl::Problem<
    1, StepperTag, FomSystemType,
    GalerkinStateType, GalerkinResidualType, GalerkinJacobianType,
    DecoderType
    >
  >
{
  using common_types = ::pressio::rom::galerkin::impl::CommonTraitsContinuousTimeApi<
    StepperTag, FomSystemType, GalerkinStateType, DecoderType>;

  using fom_system_type    = typename common_types::fom_system_type;
  using scalar_type    = typename common_types::scalar_type;
  using fom_state_type   = typename common_types::fom_state_type;
  using fom_velocity_type    = typename common_types::fom_velocity_type;
  using galerkin_state_type  = typename common_types::galerkin_state_type;
  using decoder_type   = typename common_types::decoder_type;
  using decoder_jac_type   = typename common_types::decoder_jac_type;
  using fom_state_reconstr_type  = typename common_types::fom_state_reconstr_type;
  using fom_states_manager_type  = typename common_types::fom_states_manager_type;
  using galerkin_residual_type = GalerkinResidualType;
  using galerkin_jacobian_type = GalerkinJacobianType;
  static constexpr auto binding_sentinel = common_types::binding_sentinel;
  static constexpr auto is_cont_time = true;

  using size_type = typename ::pressio::Traits<galerkin_state_type>::size_type;

  // default galerkin uses default projector, so decoderJac^T
  using projector_type = ::pressio::rom::galerkin::impl::DefaultProjector<decoder_type>;

  using residual_policy_type =
    ::pressio::rom::galerkin::impl::ResidualPolicy<
    galerkin_residual_type,
    ::pressio::rom::galerkin::impl::Projected<
      projector_type, size_type, 1,
      ::pressio::rom::galerkin::impl::DefaultFomVelocityEvaluator<
	fom_states_manager_type, fom_velocity_type, fom_system_type>
      >
    >;

  using jacobian_policy_type =
    ::pressio::rom::galerkin::impl::JacobianPolicy<
    galerkin_jacobian_type,
    ::pressio::rom::galerkin::impl::Projected<
      projector_type, size_type, 2,
      ::pressio::rom::galerkin::impl::DefaultFomApplyJacobianEvaluator<
	is_cont_time,fom_states_manager_type, decoder_jac_type, decoder_type, fom_system_type>
      >
    >;

  using stepper_type = ::pressio::ode::impl::ImplicitCompose_t<
    StepperTag, void, galerkin_state_type, residual_policy_type &, jacobian_policy_type &>;
};

// disctime implicit
template <
  std::size_t num_states,
  class FomSystemType,
  class GalerkinStateType,
  class GalerkinResidualType,
  class GalerkinJacobianType,
  class DecoderType
  >
struct Traits<
  ::pressio::rom::galerkin::impl::Problem<
    2, FomSystemType,
    GalerkinStateType, GalerkinResidualType, GalerkinJacobianType, DecoderType,
    ::pressio::ode::StepperTotalNumberOfStates<num_states>
    >
  >
{
  using common_types = ::pressio::rom::galerkin::impl::CommonTraitsDiscreteTimeApi<
    num_states, FomSystemType, GalerkinStateType, DecoderType>;

  using fom_system_type = typename common_types::fom_system_type;
  using scalar_type     = typename common_types::scalar_type;

  using fom_state_type  = typename common_types::fom_state_type;
  using fom_discrete_time_residual_type = typename common_types::fom_discrete_time_residual_type;

  using galerkin_state_type  = typename common_types::galerkin_state_type;
  using galerkin_residual_type = GalerkinResidualType;
  using galerkin_jacobian_type = GalerkinJacobianType;
  using decoder_type   = typename common_types::decoder_type;
  using decoder_jac_type   = typename common_types::decoder_jac_type;
  using fom_state_reconstr_type  = typename common_types::fom_state_reconstr_type;
  using fom_states_manager_type  = typename common_types::fom_states_manager_type;

  using size_type = typename ::pressio::Traits<galerkin_state_type>::size_type;
  static constexpr auto binding_sentinel = common_types::binding_sentinel;
  static constexpr auto is_cont_time = false;

  // default galerkin uses default projector, so decoderJac^T
  using projector_type = ::pressio::rom::galerkin::impl::DefaultProjector<decoder_type>;

  using p1 = ::pressio::rom::galerkin::impl::Projected<
    projector_type, size_type, 1,
    ::pressio::rom::galerkin::impl::FomResidualPolicyDiscreteTimeApi<
      fom_states_manager_type, fom_discrete_time_residual_type, fom_system_type
      >
    >;

  using p2 = ::pressio::rom::galerkin::impl::Projected<
    projector_type, size_type, 2,
    ::pressio::rom::galerkin::impl::DefaultFomApplyJacobianEvaluator<
      is_cont_time, fom_states_manager_type, decoder_jac_type, decoder_type, fom_system_type
      >
    >;

  using rom_system_type =
    ::pressio::rom::galerkin::impl::DiscreteTimeReducedSystem<
    scalar_type, galerkin_state_type, galerkin_residual_type, galerkin_jacobian_type,
    p1, p2
    >;

  using stepper_type = typename ::pressio::ode::impl::ImplicitComposeArb<
    num_states, const rom_system_type &, galerkin_state_type>::type;
};


//=======================
//
// MASKED
//
//=======================

// conttime explicit
template <
  class StepperTag,
  class FomSystemType,
  class GalerkinStateType,
  class DecoderType,
  class MaskerType,
  class ProjectorType
  >
struct Traits<
  ::pressio::rom::galerkin::impl::Problem<
    3, StepperTag, FomSystemType, GalerkinStateType, DecoderType, MaskerType, ProjectorType
    >
  >
{
  using common_types = ::pressio::rom::galerkin::impl::CommonTraitsContinuousTimeApi<
    StepperTag, FomSystemType, GalerkinStateType, DecoderType>;

  using fom_system_type    = typename common_types::fom_system_type;
  using scalar_type    = typename common_types::scalar_type;
  using fom_state_type   = typename common_types::fom_state_type;
  using fom_velocity_type    = typename common_types::fom_velocity_type;
  using galerkin_state_type  = typename common_types::galerkin_state_type;
  // for now, the velocity type is same as state
  using galerkin_velocity_type = galerkin_state_type;
  using decoder_type   = typename common_types::decoder_type;
  using decoder_jac_type   = typename common_types::decoder_jac_type;
  using fom_state_reconstr_type  = typename common_types::fom_state_reconstr_type;
  using fom_states_manager_type  = typename common_types::fom_states_manager_type;
  static constexpr auto binding_sentinel = common_types::binding_sentinel;
  using size_type = typename ::pressio::Traits<galerkin_state_type>::size_type;

  using masker_type = MaskerType;
  using projector_type = ProjectorType;

  static_assert
  (::pressio::rom::masker_explicit_galerkin<masker_type, scalar_type, fom_velocity_type>::value,
   "Invalid masker passed to masked Galerkin with explicit time stepping");

  static_assert
  (::pressio::rom::projector_explicit_galerkin<
   ProjectorType, scalar_type, fom_velocity_type, galerkin_velocity_type>::value,
   "Invalid projector passed to masked Galerkin with explicit time stepping");

  using rom_system_type =
    ::pressio::rom::galerkin::impl::VelocityOnlySystemUsing<
    scalar_type, galerkin_state_type, galerkin_velocity_type,
    ::pressio::rom::galerkin::impl::Projected<
      projector_type, size_type, 1,
      ::pressio::rom::galerkin::impl::Masked<
	fom_velocity_type, MaskerType,
        ::pressio::rom::galerkin::impl::DefaultFomVelocityEvaluator<
          fom_states_manager_type, fom_velocity_type, fom_system_type>
        >
      >
    >;

  // two things to fix here:
  // 1. use the public API not impl
  // 2. we need to use const rom_sytem_t & because the ROM problem
  // owns the system and the stepper only references it
  using stepper_type = typename ::pressio::ode::impl::ExplicitCompose<
    StepperTag, galerkin_state_type, const rom_system_type &>::type;
};

// conttime impiciit
template <
  class StepperTag,
  class FomSystemType,
  class GalerkinStateType,
  class GalerkinResidualType,
  class GalerkinJacobianType,
  class DecoderType,
  class MaskerType,
  class ProjectorType
  >
struct Traits<
  ::pressio::rom::galerkin::impl::Problem<
    4, StepperTag, FomSystemType,
    GalerkinStateType, GalerkinResidualType, GalerkinJacobianType,
    DecoderType, MaskerType, ProjectorType
    >
  >
{
  using common_types = ::pressio::rom::galerkin::impl::CommonTraitsContinuousTimeApi<
    StepperTag, FomSystemType, GalerkinStateType, DecoderType>;

  using fom_system_type    = typename common_types::fom_system_type;
  using scalar_type    = typename common_types::scalar_type;
  using fom_state_type   = typename common_types::fom_state_type;
  using fom_velocity_type    = typename common_types::fom_velocity_type;
  using galerkin_state_type  = typename common_types::galerkin_state_type;
  using decoder_type   = typename common_types::decoder_type;
  using decoder_jac_type   = typename common_types::decoder_jac_type;
  using fom_state_reconstr_type  = typename common_types::fom_state_reconstr_type;
  using fom_states_manager_type  = typename common_types::fom_states_manager_type;
  using galerkin_residual_type = GalerkinResidualType;
  using galerkin_jacobian_type = GalerkinJacobianType;
  static constexpr auto binding_sentinel = common_types::binding_sentinel;
  using size_type = typename ::pressio::Traits<galerkin_state_type>::size_type;
  static constexpr auto is_cont_time = true;

  using masker_type = MaskerType;
  using projector_type = ProjectorType;

  static_assert
  (::pressio::rom::masker_implicit_galerkin<
   masker_type, scalar_type, fom_velocity_type, decoder_jac_type>::value,
   "Invalid masker passed to masked Galerkin with implicit time stepping");

  static_assert
  (::pressio::rom::projector_implicit_galerkin<
   ProjectorType, scalar_type,
   fom_velocity_type, decoder_jac_type, galerkin_residual_type, galerkin_jacobian_type>::value,
   "Invalid projector passed to masked Galerkin with implicit time stepping");

  using residual_policy_type =
    ::pressio::rom::galerkin::impl::ResidualPolicy<
    galerkin_residual_type,
    ::pressio::rom::galerkin::impl::Projected<
      projector_type, size_type, 1,
      ::pressio::rom::galerkin::impl::Masked<
	fom_velocity_type, MaskerType,
        ::pressio::rom::galerkin::impl::DefaultFomVelocityEvaluator<
	  fom_states_manager_type, fom_velocity_type, fom_system_type>
        >
      >
    >;

  using jacobian_policy_type =
    ::pressio::rom::galerkin::impl::JacobianPolicy<
    galerkin_jacobian_type,
    ::pressio::rom::galerkin::impl::Projected<
      projector_type, size_type, 2,
      ::pressio::rom::galerkin::impl::Masked<
	decoder_jac_type, MaskerType,
        ::pressio::rom::galerkin::impl::DefaultFomApplyJacobianEvaluator<
	  is_cont_time, fom_states_manager_type, decoder_jac_type, decoder_type, fom_system_type>
        >
      >
    >;

  using stepper_type = ::pressio::ode::impl::ImplicitCompose_t<
    StepperTag, void, galerkin_state_type, residual_policy_type &, jacobian_policy_type &>;
};

// disctime implicit
template <
  std::size_t num_states,
  class FomSystemType,
  class GalerkinStateType,
  class GalerkinResidualType,
  class GalerkinJacobianType,
  class DecoderType,
  class ProjectorType,
  class MaskerType
  >
struct Traits<
  ::pressio::rom::galerkin::impl::Problem<
    5, FomSystemType,
    GalerkinStateType, GalerkinResidualType, GalerkinJacobianType, DecoderType, ProjectorType, MaskerType,
    ::pressio::ode::StepperTotalNumberOfStates<num_states>
    >
  >
{
  using common_types = ::pressio::rom::galerkin::impl::CommonTraitsDiscreteTimeApi<
    num_states, FomSystemType, GalerkinStateType, DecoderType>;

  using fom_system_type = typename common_types::fom_system_type;
  using scalar_type     = typename common_types::scalar_type;

  using fom_state_type  = typename common_types::fom_state_type;
  using fom_discrete_time_residual_type = typename common_types::fom_discrete_time_residual_type;

  using galerkin_state_type  = typename common_types::galerkin_state_type;
  using galerkin_residual_type = GalerkinResidualType;
  using galerkin_jacobian_type = GalerkinJacobianType;
  using decoder_type   = typename common_types::decoder_type;
  using decoder_jac_type   = typename common_types::decoder_jac_type;
  using fom_state_reconstr_type  = typename common_types::fom_state_reconstr_type;
  using fom_states_manager_type  = typename common_types::fom_states_manager_type;

  using size_type = typename ::pressio::Traits<galerkin_state_type>::size_type;
  static constexpr auto binding_sentinel = common_types::binding_sentinel;
  static constexpr auto is_cont_time = false;

  using masker_type = MaskerType;
  static_assert
  (::pressio::rom::masker_implicit_galerkin<
   masker_type, scalar_type, fom_discrete_time_residual_type, decoder_jac_type>::value,
   "Invalid masker passed to masked Galerkin with implicit time stepping");

  using projector_type = ProjectorType;
  static_assert
  (::pressio::rom::projector_implicit_galerkin<ProjectorType, scalar_type,
   fom_discrete_time_residual_type, decoder_jac_type, galerkin_residual_type, galerkin_jacobian_type>::value,
   "Invalid projector passed to hypred discrete-time Galerkin");

  using p1 = ::pressio::rom::galerkin::impl::Projected<
    projector_type, size_type, 1,
      ::pressio::rom::galerkin::impl::Masked<
	fom_discrete_time_residual_type, MaskerType,
	::pressio::rom::galerkin::impl::FomResidualPolicyDiscreteTimeApi<
	  fom_states_manager_type, fom_discrete_time_residual_type, fom_system_type
	  >
	>
    >;

  using p2 = ::pressio::rom::galerkin::impl::Projected<
    projector_type, size_type, 2,
      ::pressio::rom::galerkin::impl::Masked<
	decoder_jac_type, MaskerType,
	::pressio::rom::galerkin::impl::DefaultFomApplyJacobianEvaluator<
	  is_cont_time, fom_states_manager_type, decoder_jac_type, decoder_type, fom_system_type
	  >
	>
    >;

  using rom_system_type =
    ::pressio::rom::galerkin::impl::DiscreteTimeReducedSystem<
    scalar_type, galerkin_state_type, galerkin_residual_type, galerkin_jacobian_type,
    p1, p2>;

  using stepper_type = typename ::pressio::ode::impl::ImplicitComposeArb<
    num_states, const rom_system_type &, galerkin_state_type>::type;
};



//=======================
//
// HYP-REDUCED
//
//=======================

// contitime explicit
template <
  class StepperTag,
  class FomSystemType,
  class GalerkinStateType,
  class DecoderType,
  class ProjectorType
  >
struct Traits<
  ::pressio::rom::galerkin::impl::Problem<
    6, StepperTag, FomSystemType, GalerkinStateType, DecoderType, ProjectorType
    >
  >
{
  using common_types = ::pressio::rom::galerkin::impl::CommonTraitsContinuousTimeApi<
    StepperTag, FomSystemType, GalerkinStateType, DecoderType>;

  using fom_system_type    = typename common_types::fom_system_type;
  using scalar_type    = typename common_types::scalar_type;
  using fom_state_type   = typename common_types::fom_state_type;
  using fom_velocity_type    = typename common_types::fom_velocity_type;
  using galerkin_state_type  = typename common_types::galerkin_state_type;
  // for now, the velocity type is same as state
  using galerkin_velocity_t = galerkin_state_type;
  using decoder_type   = typename common_types::decoder_type;
  using decoder_jac_type   = typename common_types::decoder_jac_type;
  using fom_state_reconstr_type  = typename common_types::fom_state_reconstr_type;
  using fom_states_manager_type  = typename common_types::fom_states_manager_type;
  static constexpr auto binding_sentinel = common_types::binding_sentinel;
  using size_type = typename ::pressio::Traits<galerkin_state_type>::size_type;

  using projector_type = ProjectorType;
  static_assert
  (::pressio::rom::projector_explicit_galerkin<ProjectorType,
   scalar_type, fom_velocity_type, galerkin_velocity_t>::value,
   "Invalid projector passed to hypred Galerkin with explicit time stepping");

  using rom_system_type =
    ::pressio::rom::galerkin::impl::VelocityOnlySystemUsing<
    scalar_type, galerkin_state_type, galerkin_velocity_t,
    ::pressio::rom::galerkin::impl::Projected<
      projector_type, size_type, 1,
      ::pressio::rom::galerkin::impl::DefaultFomVelocityEvaluator<
        fom_states_manager_type, fom_velocity_type, fom_system_type>
      >
    >;

  // two things to fix here:
  // 1. use the public API not impl
  // 2. we need to use const rom_sytem_t & because the ROM problem
  // owns the system and the stepper only references it
  using stepper_type = typename ::pressio::ode::impl::ExplicitCompose<
    StepperTag, galerkin_state_type, const rom_system_type &>::type;
};

// conttime implicit
template <
  class StepperTag,
  class FomSystemType,
  class GalerkinStateType,
  class GalerkinResidualType,
  class GalerkinJacobianType,
  class DecoderType,
  class ProjectorType
  >
struct Traits<
  ::pressio::rom::galerkin::impl::Problem<
    7, StepperTag, FomSystemType,
    GalerkinStateType, GalerkinResidualType, GalerkinJacobianType, DecoderType, ProjectorType
    >
  >
{
  using common_types = ::pressio::rom::galerkin::impl::CommonTraitsContinuousTimeApi<
    StepperTag, FomSystemType, GalerkinStateType, DecoderType>;

  using fom_system_type    = typename common_types::fom_system_type;
  using scalar_type    = typename common_types::scalar_type;
  using fom_state_type   = typename common_types::fom_state_type;
  using fom_velocity_type    = typename common_types::fom_velocity_type;
  using galerkin_state_type  = typename common_types::galerkin_state_type;
  using decoder_type   = typename common_types::decoder_type;
  using decoder_jac_type   = typename common_types::decoder_jac_type;
  using fom_state_reconstr_type  = typename common_types::fom_state_reconstr_type;
  using fom_states_manager_type  = typename common_types::fom_states_manager_type;
  using galerkin_residual_type = GalerkinResidualType;
  using galerkin_jacobian_type = GalerkinJacobianType;
  static constexpr auto binding_sentinel = common_types::binding_sentinel;
  using size_type = typename ::pressio::Traits<galerkin_state_type>::size_type;
  static constexpr auto is_cont_time = true;

  using projector_type = ProjectorType;
  static_assert
  (::pressio::rom::projector_implicit_galerkin<ProjectorType, scalar_type,
   fom_velocity_type, decoder_jac_type, galerkin_residual_type, galerkin_jacobian_type>::value,
   "Invalid projector passed to hypred Galerkin with implicit time stepping");

  using residual_policy_type =
    ::pressio::rom::galerkin::impl::ResidualPolicy<
    galerkin_residual_type,
    ::pressio::rom::galerkin::impl::Projected<
      projector_type, size_type, 1,
      ::pressio::rom::galerkin::impl::DefaultFomVelocityEvaluator<
	fom_states_manager_type, fom_velocity_type, fom_system_type>
      >
    >;

  using jacobian_policy_type =
    ::pressio::rom::galerkin::impl::JacobianPolicy<
    galerkin_jacobian_type,
    ::pressio::rom::galerkin::impl::Projected<
      projector_type, size_type, 2,
      ::pressio::rom::galerkin::impl::DefaultFomApplyJacobianEvaluator<
	is_cont_time, fom_states_manager_type, decoder_jac_type, decoder_type, fom_system_type>
      >
    >;

  using stepper_type = ::pressio::ode::impl::ImplicitCompose_t<
    StepperTag, void, galerkin_state_type, residual_policy_type &, jacobian_policy_type &>;
};

// disctime implicit
template <
  std::size_t num_states,
  class FomSystemType,
  class GalerkinStateType,
  class GalerkinResidualType,
  class GalerkinJacobianType,
  class DecoderType,
  class ProjectorType
  >
struct Traits<
  ::pressio::rom::galerkin::impl::Problem<
    8, FomSystemType,
    GalerkinStateType, GalerkinResidualType, GalerkinJacobianType, DecoderType, ProjectorType,
    ::pressio::ode::StepperTotalNumberOfStates<num_states>
    >
  >
{
  using common_types = ::pressio::rom::galerkin::impl::CommonTraitsDiscreteTimeApi<
    num_states, FomSystemType, GalerkinStateType, DecoderType>;

  using fom_system_type = typename common_types::fom_system_type;
  using scalar_type     = typename common_types::scalar_type;

  using fom_state_type  = typename common_types::fom_state_type;
  using fom_discrete_time_residual_type = typename common_types::fom_discrete_time_residual_type;

  using galerkin_state_type  = typename common_types::galerkin_state_type;
  using galerkin_residual_type = GalerkinResidualType;
  using galerkin_jacobian_type = GalerkinJacobianType;
  using decoder_type   = typename common_types::decoder_type;
  using decoder_jac_type   = typename common_types::decoder_jac_type;
  using fom_state_reconstr_type  = typename common_types::fom_state_reconstr_type;
  using fom_states_manager_type  = typename common_types::fom_states_manager_type;

  using size_type = typename ::pressio::Traits<galerkin_state_type>::size_type;
  static constexpr auto binding_sentinel = common_types::binding_sentinel;
  static constexpr auto is_cont_time = false;

  using projector_type = ProjectorType;
  static_assert
  (::pressio::rom::projector_implicit_galerkin<ProjectorType, scalar_type,
   fom_discrete_time_residual_type, decoder_jac_type, galerkin_residual_type, galerkin_jacobian_type>::value,
   "Invalid projector passed to hypred discrete-time Galerkin");

  using p1 = ::pressio::rom::galerkin::impl::Projected<
    projector_type, size_type, 1,
    ::pressio::rom::galerkin::impl::FomResidualPolicyDiscreteTimeApi<
      fom_states_manager_type, fom_discrete_time_residual_type, fom_system_type
      >
    >;

  using p2 = ::pressio::rom::galerkin::impl::Projected<
    projector_type, size_type, 2,
    ::pressio::rom::galerkin::impl::DefaultFomApplyJacobianEvaluator<
      is_cont_time, fom_states_manager_type, decoder_jac_type, decoder_type, fom_system_type
      >
    >;

  using rom_system_type =
    ::pressio::rom::galerkin::impl::DiscreteTimeReducedSystem<
    scalar_type, galerkin_state_type, galerkin_residual_type, galerkin_jacobian_type,
    p1, p2>;

  using stepper_type = typename ::pressio::ode::impl::ImplicitComposeArb<
    num_states, const rom_system_type &, galerkin_state_type>::type;
};

}//end  namespace pressio
#endif  // ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_TRAITS_ROM_GALERKIN_COMMON_TRAITS_HPP_
