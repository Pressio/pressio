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
template <int, class ...> class ProblemContinuousTimeApi;

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
  class RomStateType,
  class decoder_type
  >
struct CommonTraitsContinuousTimeApi
{
  using fom_system_t	 = FomSystemType;
  using scalar_t	 = typename fom_system_t::scalar_type;
  using galerkin_state_t = RomStateType;

  static_assert
  (::pressio::rom::admissible_galerkin_state<galerkin_state_t>::value,
   "Invalid galerkin state type");

  using decoder_t = decoder_type;
  using decoder_jac_t = typename decoder_type::jacobian_type;

  // ---------------------
  // detect fom state type from decoder
  // ensure it is consistent with the fom_state_type from the app
  using fom_state_from_decoder_t = typename decoder_type::fom_state_type;
  using fom_state_from_adapter_t = typename fom_system_t::state_type;
  static_assert
  (std::is_same<fom_state_from_decoder_t, fom_state_from_adapter_t>::value,
   "The fom state type detected from the fom adapter must match the fom state type used in the decoder");
  using fom_state_t = fom_state_from_decoder_t;

  // ---------------------
  // we don't allow fom state and fom velocity to have different types
  // but need to make sure this assumption is consistent with fom class
  using fom_velocity_t = fom_state_t;
  using fom_velocity_from_adapter_t = typename fom_system_t::velocity_type;
  static_assert
  (std::is_same<fom_state_t, fom_velocity_from_adapter_t>::value,
   "Currently, the fom state and velocity must be of the same type");

  // ---------------------
  // fom state reconstructor
  using fom_state_reconstr_t = ::pressio::rom::FomStateReconstructor<decoder_t>;

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
  using fom_states_manager_t = mpl::conditional_t<
    ::pressio::ode::is_explicit_stepper_tag<OdeTag>::value,
    ::pressio::rom::ManagerFomStatesUnsteadyExplicit<
      fom_state_t, fom_state_reconstr_t, nstates>,
    ::pressio::rom::ManagerFomStatesUnsteadyImplicit<
      fom_state_t, fom_state_reconstr_t, nstates>
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

}}} // end namespace pressio::rom::galekin::impl

//
// DEFAULT problem
//
template <
  class StepperTag,
  class FomSystemType,
  class RomStateType,
  class DecoderType
  >
struct Traits<
  ::pressio::rom::galerkin::impl::ProblemContinuousTimeApi<
    0, StepperTag, FomSystemType, RomStateType, DecoderType
    >
  >
{
  using common_types_t = ::pressio::rom::galerkin::impl::CommonTraitsContinuousTimeApi<
    StepperTag, FomSystemType, RomStateType, DecoderType>;
  using fom_system_t		= typename common_types_t::fom_system_t;
  using scalar_t		= typename common_types_t::scalar_t;
  using fom_state_t		= typename common_types_t::fom_state_t;
  using fom_velocity_t		= typename common_types_t::fom_velocity_t;
  using galerkin_state_t	= typename common_types_t::galerkin_state_t;
  // for now, the velocity type is same as state
  using galerkin_velocity_t = galerkin_state_t;

  using decoder_t		= typename common_types_t::decoder_t;
  using decoder_jac_t		= typename common_types_t::decoder_jac_t;
  using fom_state_reconstr_t	= typename common_types_t::fom_state_reconstr_t;
  using fom_states_manager_t	= typename common_types_t::fom_states_manager_t;
  static constexpr auto binding_sentinel = common_types_t::binding_sentinel;
  using size_type = typename ::pressio::Traits<galerkin_state_t>::size_type;

  // default galerkin uses default projector, so decoderJac^T
  using projector_t = ::pressio::rom::galerkin::impl::DefaultProjector<decoder_t>;

  using rom_system_t =
    ::pressio::rom::galerkin::impl::VelocityOnlySystemUsing<
    scalar_t, galerkin_state_t, galerkin_velocity_t,
    ::pressio::rom::galerkin::impl::Projected<
      projector_t, size_type, 1,
      ::pressio::rom::galerkin::impl::DefaultFomVelocityEvaluator<
        fom_states_manager_t, fom_velocity_t, fom_system_t>
      >
    >;

  // two things to fix here:
  // 1. use the public API not impl
  // 2. we need to use const rom_sytem_t & because the ROM problem
  // owns the system and the stepper only references it
  using stepper_t = typename ::pressio::ode::impl::ExplicitCompose<
    StepperTag, galerkin_state_t, const rom_system_t &>::type;
};

template <
  class StepperTag,
  class FomSystemType,
  class RomStateType,
  class RomResidualType,
  class RomJacobianType,
  class DecoderType
  >
struct Traits<
  ::pressio::rom::galerkin::impl::ProblemContinuousTimeApi<
    1, StepperTag, FomSystemType,
    RomStateType, RomResidualType, RomJacobianType,
    DecoderType
    >
  >
{
  using common_types_t = ::pressio::rom::galerkin::impl::CommonTraitsContinuousTimeApi<
    StepperTag, FomSystemType, RomStateType, DecoderType>;

  using fom_system_t    = typename common_types_t::fom_system_t;
  using scalar_t    = typename common_types_t::scalar_t;
  using fom_state_t   = typename common_types_t::fom_state_t;
  using fom_velocity_t    = typename common_types_t::fom_velocity_t;
  using galerkin_state_t  = typename common_types_t::galerkin_state_t;
  using decoder_t   = typename common_types_t::decoder_t;
  using decoder_jac_t   = typename common_types_t::decoder_jac_t;
  using fom_state_reconstr_t  = typename common_types_t::fom_state_reconstr_t;
  using fom_states_manager_t  = typename common_types_t::fom_states_manager_t;
  using galerkin_residual_t = RomResidualType;
  using galerkin_jacobian_t = RomJacobianType;
  static constexpr auto binding_sentinel = common_types_t::binding_sentinel;

  using size_type = typename ::pressio::Traits<galerkin_state_t>::size_type;

  // default galerkin uses default projector, so decoderJac^T
  using projector_t = ::pressio::rom::galerkin::impl::DefaultProjector<decoder_t>;

  using residual_policy_t =
    ::pressio::rom::galerkin::impl::ResidualPolicy<
    galerkin_residual_t,
    ::pressio::rom::galerkin::impl::Projected<
      projector_t, size_type, 1,
      ::pressio::rom::galerkin::impl::DefaultFomVelocityEvaluator<
	fom_states_manager_t, fom_velocity_t, fom_system_t>
      >
    >;

  using jacobian_policy_t =
    ::pressio::rom::galerkin::impl::JacobianPolicy<
    galerkin_jacobian_t,
    ::pressio::rom::galerkin::impl::Projected<
      projector_t, size_type, 2,
      ::pressio::rom::galerkin::impl::DefaultFomApplyJacobianEvaluator<
	fom_states_manager_t, decoder_jac_t, decoder_t, fom_system_t>
      >
    >;

  using stepper_t = ::pressio::ode::impl::ImplicitCompose_t<
    StepperTag, void, galerkin_state_t, residual_policy_t &, jacobian_policy_t &>;
};


//
// MASKED VELO problem
//
template <
  class StepperTag,
  class FomSystemType,
  class RomStateType,
  class DecoderType,
  class MaskerType,
  class ProjectorType
  >
struct Traits<
  ::pressio::rom::galerkin::impl::ProblemContinuousTimeApi<
    2, StepperTag, FomSystemType, RomStateType, DecoderType, MaskerType, ProjectorType
    >
  >
{
  using common_types_t = ::pressio::rom::galerkin::impl::CommonTraitsContinuousTimeApi<
    StepperTag, FomSystemType, RomStateType, DecoderType>;

  using fom_system_t    = typename common_types_t::fom_system_t;
  using scalar_t    = typename common_types_t::scalar_t;
  using fom_state_t   = typename common_types_t::fom_state_t;
  using fom_velocity_t    = typename common_types_t::fom_velocity_t;
  using galerkin_state_t  = typename common_types_t::galerkin_state_t;
  // for now, the velocity type is same as state
  using galerkin_velocity_t = galerkin_state_t;
  using decoder_t   = typename common_types_t::decoder_t;
  using decoder_jac_t   = typename common_types_t::decoder_jac_t;
  using fom_state_reconstr_t  = typename common_types_t::fom_state_reconstr_t;
  using fom_states_manager_t  = typename common_types_t::fom_states_manager_t;
  static constexpr auto binding_sentinel = common_types_t::binding_sentinel;
  using size_type = typename ::pressio::Traits<galerkin_state_t>::size_type;

  using masker_t = MaskerType;
  using projector_t = ProjectorType;

  static_assert
  (::pressio::rom::masker_explicit_galerkin<masker_t, scalar_t, fom_velocity_t>::value,
   "Invalid masker passed to masked Galerkin with explicit time stepping");

  static_assert
  (::pressio::rom::projector_explicit_galerkin<ProjectorType,
      fom_velocity_t, galerkin_velocity_t>::value,
   "Invalid projector passed to masked Galerkin with explicit time stepping");

  using rom_system_t =
    ::pressio::rom::galerkin::impl::VelocityOnlySystemUsing<
    scalar_t, galerkin_state_t, galerkin_velocity_t,
    ::pressio::rom::galerkin::impl::Projected<
      projector_t, size_type, 1,
      ::pressio::rom::galerkin::impl::Masked<
	fom_velocity_t, MaskerType,
        ::pressio::rom::galerkin::impl::DefaultFomVelocityEvaluator<
          fom_states_manager_t, fom_velocity_t, fom_system_t>
        >
      >
    >;

  // two things to fix here:
  // 1. use the public API not impl
  // 2. we need to use const rom_sytem_t & because the ROM problem
  // owns the system and the stepper only references it
  using stepper_t = typename ::pressio::ode::impl::ExplicitCompose<
    StepperTag, galerkin_state_t, const rom_system_t &>::type;
};

template <
  class StepperTag,
  class FomSystemType,
  class RomStateType,
  class RomResidualType,
  class RomJacobianType,
  class DecoderType,
  class MaskerType,
  class ProjectorType
  >
struct Traits<
  ::pressio::rom::galerkin::impl::ProblemContinuousTimeApi<
    3, StepperTag, FomSystemType,
    RomStateType, RomResidualType, RomJacobianType,
    DecoderType, MaskerType, ProjectorType
    >
  >
{
  using common_types_t = ::pressio::rom::galerkin::impl::CommonTraitsContinuousTimeApi<
    StepperTag, FomSystemType, RomStateType, DecoderType>;

  using fom_system_t    = typename common_types_t::fom_system_t;
  using scalar_t    = typename common_types_t::scalar_t;
  using fom_state_t   = typename common_types_t::fom_state_t;
  using fom_velocity_t    = typename common_types_t::fom_velocity_t;
  using galerkin_state_t  = typename common_types_t::galerkin_state_t;
  using decoder_t   = typename common_types_t::decoder_t;
  using decoder_jac_t   = typename common_types_t::decoder_jac_t;
  using fom_state_reconstr_t  = typename common_types_t::fom_state_reconstr_t;
  using fom_states_manager_t  = typename common_types_t::fom_states_manager_t;
  using galerkin_residual_t = RomResidualType;
  using galerkin_jacobian_t = RomJacobianType;
  static constexpr auto binding_sentinel = common_types_t::binding_sentinel;
  using size_type = typename ::pressio::Traits<galerkin_state_t>::size_type;

  using masker_t = MaskerType;
  using projector_t = ProjectorType;

  using residual_policy_t =
    ::pressio::rom::galerkin::impl::ResidualPolicy<
    galerkin_residual_t,
    ::pressio::rom::galerkin::impl::Projected<
      projector_t, size_type, 1,
      ::pressio::rom::galerkin::impl::Masked<
	fom_velocity_t, MaskerType,
        ::pressio::rom::galerkin::impl::DefaultFomVelocityEvaluator<
	  fom_states_manager_t, fom_velocity_t, fom_system_t>
        >
      >
    >;

  using jacobian_policy_t =
    ::pressio::rom::galerkin::impl::JacobianPolicy<
    galerkin_jacobian_t,
    ::pressio::rom::galerkin::impl::Projected<
      projector_t, size_type, 2,
      ::pressio::rom::galerkin::impl::Masked<
	decoder_jac_t, MaskerType,
        ::pressio::rom::galerkin::impl::DefaultFomApplyJacobianEvaluator<
	  fom_states_manager_t, decoder_jac_t, decoder_t, fom_system_t>
        >
      >
    >;

  using stepper_t = ::pressio::ode::impl::ImplicitCompose_t<
    StepperTag, void, galerkin_state_t, residual_policy_t &, jacobian_policy_t &>;
};


//
// HYPRED VELO problem
//
template <
  class StepperTag,
  class FomSystemType,
  class RomStateType,
  class DecoderType,
  class ProjectorType
  >
struct Traits<
  ::pressio::rom::galerkin::impl::ProblemContinuousTimeApi<
    4, StepperTag, FomSystemType, RomStateType, DecoderType, ProjectorType
    >
  >
{
  using common_types_t = ::pressio::rom::galerkin::impl::CommonTraitsContinuousTimeApi<
    StepperTag, FomSystemType, RomStateType, DecoderType>;

  using fom_system_t    = typename common_types_t::fom_system_t;
  using scalar_t    = typename common_types_t::scalar_t;
  using fom_state_t   = typename common_types_t::fom_state_t;
  using fom_velocity_t    = typename common_types_t::fom_velocity_t;
  using galerkin_state_t  = typename common_types_t::galerkin_state_t;
  // for now, the velocity type is same as state
  using galerkin_velocity_t = galerkin_state_t;
  using decoder_t   = typename common_types_t::decoder_t;
  using decoder_jac_t   = typename common_types_t::decoder_jac_t;
  using fom_state_reconstr_t  = typename common_types_t::fom_state_reconstr_t;
  using fom_states_manager_t  = typename common_types_t::fom_states_manager_t;
  static constexpr auto binding_sentinel = common_types_t::binding_sentinel;
  using size_type = typename ::pressio::Traits<galerkin_state_t>::size_type;

  using projector_t = ProjectorType;
  static_assert
  (::pressio::rom::projector_explicit_galerkin<ProjectorType,
      fom_velocity_t, galerkin_velocity_t>::value,
   "Invalid projector passed to hypred velo Galerkin with explicit time stepping");

  using rom_system_t =
    ::pressio::rom::galerkin::impl::VelocityOnlySystemUsing<
    scalar_t, galerkin_state_t, galerkin_velocity_t,
    ::pressio::rom::galerkin::impl::Projected<
      projector_t, size_type, 1,
      ::pressio::rom::galerkin::impl::DefaultFomVelocityEvaluator<
        fom_states_manager_t, fom_velocity_t, fom_system_t>
      >
    >;

  // two things to fix here:
  // 1. use the public API not impl
  // 2. we need to use const rom_sytem_t & because the ROM problem
  // owns the system and the stepper only references it
  using stepper_t = typename ::pressio::ode::impl::ExplicitCompose<
    StepperTag, galerkin_state_t, const rom_system_t &>::type;
};

template <
  class StepperTag,
  class FomSystemType,
  class RomStateType,
  class RomResidualType,
  class RomJacobianType,
  class DecoderType,
  class ProjectorType
  >
struct Traits<
  ::pressio::rom::galerkin::impl::ProblemContinuousTimeApi<
    5, StepperTag, FomSystemType,
    RomStateType, RomResidualType, RomJacobianType, DecoderType, ProjectorType
    >
  >
{
  using common_types_t = ::pressio::rom::galerkin::impl::CommonTraitsContinuousTimeApi<
    StepperTag, FomSystemType, RomStateType, DecoderType>;

  using fom_system_t    = typename common_types_t::fom_system_t;
  using scalar_t    = typename common_types_t::scalar_t;
  using fom_state_t   = typename common_types_t::fom_state_t;
  using fom_velocity_t    = typename common_types_t::fom_velocity_t;
  using galerkin_state_t  = typename common_types_t::galerkin_state_t;
  using decoder_t   = typename common_types_t::decoder_t;
  using decoder_jac_t   = typename common_types_t::decoder_jac_t;
  using fom_state_reconstr_t  = typename common_types_t::fom_state_reconstr_t;
  using fom_states_manager_t  = typename common_types_t::fom_states_manager_t;
  using galerkin_residual_t = RomResidualType;
  using galerkin_jacobian_t = RomJacobianType;
  static constexpr auto binding_sentinel = common_types_t::binding_sentinel;
  using size_type = typename ::pressio::Traits<galerkin_state_t>::size_type;

  using projector_t = ProjectorType;

  using residual_policy_t =
    ::pressio::rom::galerkin::impl::ResidualPolicy<
    galerkin_residual_t,
    ::pressio::rom::galerkin::impl::Projected<
      projector_t, size_type, 1,
      ::pressio::rom::galerkin::impl::DefaultFomVelocityEvaluator<
	fom_states_manager_t, fom_velocity_t, fom_system_t>
      >
    >;

  using jacobian_policy_t =
    ::pressio::rom::galerkin::impl::JacobianPolicy<
    galerkin_jacobian_t,
    ::pressio::rom::galerkin::impl::Projected<
      projector_t, size_type, 2,
      ::pressio::rom::galerkin::impl::DefaultFomApplyJacobianEvaluator<
	fom_states_manager_t, decoder_jac_t, decoder_t, fom_system_t>
      >
    >;

  using stepper_t = ::pressio::ode::impl::ImplicitCompose_t<
    StepperTag, void, galerkin_state_t, residual_policy_t &, jacobian_policy_t &>;
};


}//end  namespace pressio
#endif  // ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_TRAITS_ROM_GALERKIN_COMMON_TRAITS_HPP_
