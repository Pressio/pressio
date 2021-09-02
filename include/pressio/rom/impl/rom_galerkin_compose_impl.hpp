/*
//@HEADER
// ************************************************************************
//
// rom_compose_impl.hpp
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

#ifndef ROM_GALERKIN_IMPL_ROM_COMPOSE_IMPL_HPP_
#define ROM_GALERKIN_IMPL_ROM_COMPOSE_IMPL_HPP_

#include "./rom_problem_members_common_mixins.hpp"
#include "./rom_galerkin_problem_members.hpp"
#include "./rom_galerkin_decorator_projected.hpp"
#include "./rom_galerkin_decorator_masked.hpp"
#include "./rom_galerkin_fom_velocity_policy.hpp"
#include "./rom_galerkin_fom_residual_policy.hpp"
#include "./rom_galerkin_fom_apply_jacobian_policy.hpp"
#include "./rom_galerkin_systems.hpp"
#include "./rom_galerkin_policy_residual.hpp"
#include "./rom_galerkin_policy_jacobian.hpp"
#include "./rom_galerkin_problem_traits.hpp"
#include "./rom_galerkin_problem.hpp"

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template<typename T, class = void>
struct select_galerkin_types{
  using residual_type = void;
  using jacobian_type = void;
};

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<typename T>
struct select_galerkin_types<
  T, mpl::enable_if_t<::pressio::is_dynamic_vector_eigen<T>::value>
  >
{
  // for now use residual_type = state_type
  using residual_type = T;
  using jacobian_type = Eigen::Matrix<typename Traits<T>::scalar_type, -1,-1>;
};
#endif

template <typename tag>
struct supported_explicit_stepper_tag{
  static_assert
  (std::is_same<tag, ::pressio::ode::ForwardEuler>::value or
   std::is_same<tag, ::pressio::ode::RungeKutta4>::value or
   std::is_same<tag, ::pressio::ode::SSPRungeKutta3>::value or
   std::is_same<tag, ::pressio::ode::AdamsBashforth2>::value,
   "The explicit stepper tag you are passing to create the galerkin problem \
is not supported: this can be because the Galerkin implementation does \
not support it, or because you added a new ode scheme in the ode package \
but forgot to update the list of explicit tags supported by Galerkin which");

  static constexpr auto value = true;
};

template <typename tag>
struct supported_implicit_stepper_tag{
  static_assert
  (std::is_same<tag, ::pressio::ode::BDF1>::value or
   std::is_same<tag, ::pressio::ode::BDF2>::value or
   std::is_same<tag, ::pressio::ode::CrankNicolson>::value,
   "The implicit stepper tag you are passing to create the galerkin problem \
is not supported: this can be because the Galerkin implementation does \
not support it, or because you added a new ode scheme in the ode package \
but forgot to update the list of implicit tags supported by Galerkin");

  static constexpr auto value = true;
};


/*
  NOTE (A): pay attention below since for implicit time stepping we need to
  know the type of the reduced residual and jacobian.
  Obviously, this should be something that is compatible
  with the galerkin state and such that we can solve the reduced problem
  easily: typically, the galerkin state is a shared-mem array, e.g. eigen or Kokkos.
  So when the user does not pass them explicitly, we check what type is
  the galerkin_state_t and from that set the types for the galerkin residual and jacobian.
  For example, if galerkin_state = eigen_vector, then it would make sense
  to set galerkin_residual = galerkin_state and galerkin_jacobian = eigen_matrix
*/

template<class stepper_tag, class ...Args>
struct ComposeDefaultProblemContTime;

template<class stepper_tag, class ...Args>
struct ComposeHypRedVeloProblemContTime;

template<class stepper_tag, class ...Args>
struct ComposeMaskedProblemContTime;

//////////////////////////////////////
/// default
//////////////////////////////////////

/// explicit stepping
template<
  class StepperTag,
  class FomSystemType,
  class DecoderType,
  class GalerkinStateType,
  class FomReferenceState
  >
struct ComposeDefaultProblemContTime<
  StepperTag,
  mpl::enable_if_t< ::pressio::ode::is_explicit_stepper_tag<StepperTag>::value >,
  FomSystemType, DecoderType, GalerkinStateType, FomReferenceState>
{
  static_assert(supported_explicit_stepper_tag<StepperTag>::value, "");

  static_assert
  (::pressio::rom::continuous_time_fom_system_with_at_least_velocity<
   mpl::remove_cvref_t<FomSystemType>>::value,
   "The FOM system does not satisfy the continuous-time concept with velocity");

  static_assert(::pressio::rom::decoder<mpl::remove_cvref_t<DecoderType>, GalerkinStateType>::value,
		"Invalid decoder detected");

  static_assert
    (std::is_same<mpl::remove_cvref_t<FomReferenceState>,
     typename mpl::remove_cvref_t<DecoderType>::fom_state_type>::value,
   "The type deduced for the FOM nominal state passed to the create function is not \
   compatible with the FOM state type detected from adapter class");

  using type = ::pressio::rom::galerkin::impl::Problem<
    0, StepperTag, FomSystemType, GalerkinStateType, DecoderType>;
};

/// implicit stepping
template<
  class StepperTag,
  class FomSystemType,
  class DecoderType,
  class GalerkinStateType,
  class FomReferenceState
  >
struct ComposeDefaultProblemContTime<
  StepperTag,
  mpl::enable_if_t<::pressio::ode::is_implicit_stepper_tag<StepperTag>::value >,
  FomSystemType, DecoderType, GalerkinStateType, FomReferenceState>
{
  static_assert(supported_implicit_stepper_tag<StepperTag>::value, "");

  static_assert(::pressio::rom::decoder<mpl::remove_cvref_t<DecoderType>, GalerkinStateType>::value,
    "Invalid decoder detected");
  using decoder_jacobian_type = typename mpl::remove_cvref_t<DecoderType>::jacobian_type;

  static_assert
  (::pressio::rom::continuous_time_fom_system_with_user_provided_apply_jacobian<
   mpl::remove_cvref_t<FomSystemType>, decoder_jacobian_type>::value,
   "The FOM system does not satisfy the continuous-time concept with velocity and apply jacobian");

  static_assert
    (std::is_same<mpl::remove_cvref_t<FomReferenceState>,
     typename mpl::remove_cvref_t<DecoderType>::fom_state_type>::value,
   "The type deduced for the FOM nominal state passed to the create function is not \
   compatible with the FOM state type detected from decoder");

  // infer residual and jacobian from state, see NOTE at top
  using galerkin_residual_t = typename select_galerkin_types<GalerkinStateType>::residual_type;
  using galerkin_jacobian_t = typename select_galerkin_types<GalerkinStateType>::jacobian_type;

  using type = ::pressio::rom::galerkin::impl::Problem<
    1, StepperTag, FomSystemType, GalerkinStateType,
    galerkin_residual_t, galerkin_jacobian_t, DecoderType>;
};

template<
  std::size_t num_states,
  class FomSystemType,
  class DecoderType,
  class GalerkinStateType,
  class FomReferenceState
  >
struct ComposeDefaultProblemDiscTime
{
  static_assert(::pressio::rom::decoder<mpl::remove_cvref_t<DecoderType>, GalerkinStateType>::value,
    "Invalid decoder detected");
  using decoder_jacobian_type = typename mpl::remove_cvref_t<DecoderType>::jacobian_type;

  static_assert
  (::pressio::rom::discrete_time_fom_system_with_user_provided_apply_jacobian<
   mpl::remove_cvref_t<FomSystemType>, num_states, decoder_jacobian_type>::value,
   "Discrete-time Galerkin: FOM system must meet the discrete-time API with jacobian");

  static_assert
    (std::is_same<mpl::remove_cvref_t<FomReferenceState>,
     typename mpl::remove_cvref_t<DecoderType>::fom_state_type>::value,
   "The type deduced for the FOM nominal state passed to the create function is not \
   compatible with the FOM state type detected from decoder");

  // infer residual and jacobian from state, see NOTE at top
  using galerkin_residual_t = typename select_galerkin_types<GalerkinStateType>::residual_type;
  using galerkin_jacobian_t = typename select_galerkin_types<GalerkinStateType>::jacobian_type;

  using type = ::pressio::rom::galerkin::impl::Problem<
    2, FomSystemType, GalerkinStateType,
    galerkin_residual_t, galerkin_jacobian_t, DecoderType,
    ::pressio::ode::StepperTotalNumberOfStates<num_states>>;
};


//////////////////////////////////////
/// masked
//////////////////////////////////////

// explicit stepping
template<
  class StepperTag,
  class FomSystemType,
  class DecoderType,
  class GalerkinStateType,
  class FomReferenceState,
  class MaskerType,
  class ProjectorType
  >
struct ComposeMaskedProblemContTime<
  StepperTag,
  mpl::enable_if_t< ::pressio::ode::is_explicit_stepper_tag<StepperTag>::value >,
  FomSystemType, DecoderType, GalerkinStateType, FomReferenceState, MaskerType, ProjectorType>
{
  static_assert(supported_explicit_stepper_tag<StepperTag>::value, "");

  static_assert
  (::pressio::rom::continuous_time_fom_system_with_at_least_velocity<
   mpl::remove_cvref_t<FomSystemType>>::value,
   "The FOM system does not satisfy the continuous-time concept with velocity");

  static_assert(::pressio::rom::decoder<mpl::remove_cvref_t<DecoderType>, GalerkinStateType>::value,
    "Invalid decoder detected");

  static_assert
    (std::is_same<mpl::remove_cvref_t<FomReferenceState>,
     typename mpl::remove_cvref_t<DecoderType>::fom_state_type>::value,
   "The type deduced for the FOM nominal state passed to the create function is not \
   compatible with the FOM state type detected from decoder");

  using type = ::pressio::rom::galerkin::impl::Problem<
    3, StepperTag, FomSystemType, GalerkinStateType, DecoderType, MaskerType, ProjectorType>;
};

/// implicit stepping
template<
  class StepperTag,
  class FomSystemType,
  class DecoderType,
  class GalerkinStateType,
  class FomReferenceState,
  class MaskerType,
  class ProjectorType
  >
struct ComposeMaskedProblemContTime<
  StepperTag,
  mpl::enable_if_t<::pressio::ode::is_implicit_stepper_tag<StepperTag>::value >,
  FomSystemType, DecoderType, GalerkinStateType, FomReferenceState, MaskerType, ProjectorType>
{
  static_assert(supported_implicit_stepper_tag<StepperTag>::value, "");

  static_assert(::pressio::rom::decoder<mpl::remove_cvref_t<DecoderType>, GalerkinStateType>::value,
    "Invalid decoder detected");
  using decoder_jacobian_type = typename mpl::remove_cvref_t<DecoderType>::jacobian_type;

  static_assert
  (::pressio::rom::continuous_time_fom_system_with_user_provided_apply_jacobian<
   mpl::remove_cvref_t<FomSystemType>, decoder_jacobian_type>::value,
   "The FOM system does not satisfy the continuous-time concept with velocity and apply jacobian");

  static_assert
    (std::is_same<mpl::remove_cvref_t<FomReferenceState>,
     typename mpl::remove_cvref_t<DecoderType>::fom_state_type>::value,
   "The type deduced for the FOM nominal state passed to the create function is not \
   compatible with the FOM state type detected from decoder");

  // infer residual and jacobian from state, see NOTE at top
  using galerkin_residual_t = typename select_galerkin_types<GalerkinStateType>::residual_type;
  using galerkin_jacobian_t = typename select_galerkin_types<GalerkinStateType>::jacobian_type;

  using type = ::pressio::rom::galerkin::impl::Problem<
    4, StepperTag, FomSystemType, GalerkinStateType,
    galerkin_residual_t, galerkin_jacobian_t, DecoderType,
    MaskerType, ProjectorType>;
};

template<
  std::size_t num_states,
  class FomSystemType,
  class DecoderType,
  class GalerkinStateType,
  class FomReferenceState,
  class ProjectorType,
  class MaskerType
  >
struct ComposeMaskedProblemDiscTime
{
  static_assert
  (::pressio::rom::decoder<mpl::remove_cvref_t<DecoderType>, GalerkinStateType>::value,
   "Invalid decoder detected");
  using decoder_jacobian_type = typename mpl::remove_cvref_t<DecoderType>::jacobian_type;

  static_assert
  (::pressio::rom::discrete_time_fom_system_with_user_provided_apply_jacobian<
   mpl::remove_cvref_t<FomSystemType>, num_states, decoder_jacobian_type>::value,
   "Discrete-time Galerkin: FOM system must meet the discrete-time API with jacobian");

  static_assert
    (std::is_same<mpl::remove_cvref_t<FomReferenceState>,
     typename mpl::remove_cvref_t<DecoderType>::fom_state_type>::value,
   "The type deduced for the FOM nominal state passed to the create function is not \
   compatible with the FOM state type detected from decoder");

  // infer residual and jacobian from state, see NOTE at top
  using galerkin_residual_t = typename select_galerkin_types<GalerkinStateType>::residual_type;
  using galerkin_jacobian_t = typename select_galerkin_types<GalerkinStateType>::jacobian_type;

  using type = ::pressio::rom::galerkin::impl::Problem<
    5, FomSystemType, GalerkinStateType,
    galerkin_residual_t, galerkin_jacobian_t, DecoderType, ProjectorType, MaskerType,
    ::pressio::ode::StepperTotalNumberOfStates<num_states>>;
};


//////////////////////////////////////
/// hyper-reduced
//////////////////////////////////////

// explicit stepping
template<
  class StepperTag,
  class FomSystemType,
  class DecoderType,
  class GalerkinStateType,
  class FomReferenceState,
  class ProjectorType
  >
struct ComposeHypRedVeloProblemContTime<
  StepperTag,
  mpl::enable_if_t< ::pressio::ode::is_explicit_stepper_tag<StepperTag>::value >,
  FomSystemType, DecoderType, GalerkinStateType, FomReferenceState, ProjectorType>
{
  static_assert(supported_explicit_stepper_tag<StepperTag>::value, "");

  static_assert
  (::pressio::rom::continuous_time_fom_system_with_at_least_velocity<
   mpl::remove_cvref_t<FomSystemType>>::value,
   "The FOM system does not satisfy the continuous-time concept with velocity");

  static_assert(::pressio::rom::decoder<mpl::remove_cvref_t<DecoderType>, GalerkinStateType>::value,
    "Invalid decoder detected");

  static_assert
    (std::is_same<mpl::remove_cvref_t<FomReferenceState>,
     typename mpl::remove_cvref_t<DecoderType>::fom_state_type>::value,
   "The type deduced for the FOM nominal state passed to the create function is not \
   compatible with the FOM state type detected from decoder");

  using type = ::pressio::rom::galerkin::impl::Problem<
    6, StepperTag, FomSystemType, GalerkinStateType, DecoderType, ProjectorType>;
};

/// implicit stepping
template<
  class StepperTag,
  class FomSystemType,
  class DecoderType,
  class GalerkinStateType,
  class FomReferenceState,
  class ProjectorType
  >
struct ComposeHypRedVeloProblemContTime<
  StepperTag,
  mpl::enable_if_t<::pressio::ode::is_implicit_stepper_tag<StepperTag>::value >,
  FomSystemType, DecoderType, GalerkinStateType, FomReferenceState, ProjectorType>
{
  static_assert(supported_implicit_stepper_tag<StepperTag>::value, "");

  static_assert(::pressio::rom::decoder<mpl::remove_cvref_t<DecoderType>, GalerkinStateType>::value,
    "Invalid decoder detected");
  using decoder_jacobian_type = typename mpl::remove_cvref_t<DecoderType>::jacobian_type;

  static_assert
  (::pressio::rom::continuous_time_fom_system_with_user_provided_apply_jacobian<
   mpl::remove_cvref_t<FomSystemType>, decoder_jacobian_type>::value,
   "The FOM system does not satisfy the continuous-time concept with velocity and apply jacobian");

  static_assert
    (std::is_same<mpl::remove_cvref_t<FomReferenceState>,
     typename mpl::remove_cvref_t<DecoderType>::fom_state_type>::value,
   "The type deduced for the FOM nominal state passed to the create function is not \
   compatible with the FOM state type detected from decoder");

  // infer residual and jacobian from state, see NOTE at top
  using galerkin_residual_t = typename select_galerkin_types<GalerkinStateType>::residual_type;
  using galerkin_jacobian_t = typename select_galerkin_types<GalerkinStateType>::jacobian_type;

  using type = ::pressio::rom::galerkin::impl::Problem<
    7, StepperTag, FomSystemType, GalerkinStateType,
    galerkin_residual_t, galerkin_jacobian_t, DecoderType, ProjectorType>;
};

template<
  std::size_t num_states,
  class FomSystemType,
  class DecoderType,
  class GalerkinStateType,
  class FomReferenceState,
  class ProjectorType

  >
struct ComposeHypRedProblemDiscTime
{
  static_assert(::pressio::rom::decoder<mpl::remove_cvref_t<DecoderType>, GalerkinStateType>::value,
    "Invalid decoder detected");
  using decoder_jacobian_type = typename mpl::remove_cvref_t<DecoderType>::jacobian_type;

  static_assert
  (::pressio::rom::discrete_time_fom_system_with_user_provided_apply_jacobian<
   mpl::remove_cvref_t<FomSystemType>, num_states, decoder_jacobian_type>::value,
   "Discrete-time Galerkin: FOM system must meet the discrete-time API with jacobian");

  static_assert
    (std::is_same<mpl::remove_cvref_t<FomReferenceState>,
     typename mpl::remove_cvref_t<DecoderType>::fom_state_type>::value,
   "The type deduced for the FOM nominal state passed to the create function is not \
   compatible with the FOM state type detected from decoder");

  // infer residual and jacobian from state, see NOTE at top
  using galerkin_residual_t = typename select_galerkin_types<GalerkinStateType>::residual_type;
  using galerkin_jacobian_t = typename select_galerkin_types<GalerkinStateType>::jacobian_type;

  using type = ::pressio::rom::galerkin::impl::Problem<
    8, FomSystemType, GalerkinStateType,
    galerkin_residual_t, galerkin_jacobian_t, DecoderType, ProjectorType,
    ::pressio::ode::StepperTotalNumberOfStates<num_states>>;
};


template<typename ...Args>
using ComposeDefaultProblemContTime_t =
  typename ComposeDefaultProblemContTime<Args...>::type;

template<typename ...Args>
using ComposeMaskedVelocityProblemContTime_t =
  typename ComposeMaskedProblemContTime<Args...>::type;

template<typename ...Args>
using ComposeHypRedVeloProblemContTime_t =
  typename ComposeHypRedVeloProblemContTime<Args...>::type;

}}}}
#endif  // ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_ROM_COMPOSE_IMPL_HPP_
