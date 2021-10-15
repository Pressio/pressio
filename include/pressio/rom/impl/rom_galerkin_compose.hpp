/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_compose.hpp
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

#ifndef ROM_IMPL_ROM_GALERKIN_COMPOSE_HPP_
#define ROM_IMPL_ROM_GALERKIN_COMPOSE_HPP_

#include "./rom_problem_members_common_mixins.hpp"
#include "./rom_galerkin_problem_members.hpp"
#include "./rom_galerkin_projector_default.hpp"
#include "./rom_galerkin_decorator_projected.hpp"
#include "./rom_galerkin_decorator_masked.hpp"
#include "./rom_galerkin_fom_velocity_policy.hpp"
#include "./rom_galerkin_fom_residual_policy.hpp"
#include "./rom_galerkin_fom_apply_jacobian_policy.hpp"
#include "./rom_galerkin_systems.hpp"
#include "./rom_galerkin_policy_residual.hpp"
#include "./rom_galerkin_policy_jacobian.hpp"
#include "./rom_galerkin_traits.hpp"
#include "./rom_galerkin_problem.hpp"

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

/*
  NOTE: pay attention below since for implicit time stepping we need to
  know the type of the reduced residual and jacobian.
  Obviously, this should be something that is compatible
  with the galerkin state and such that we can solve the reduced problem
  easily: typically, the galerkin state is a shared-mem array, e.g. eigen or Kokkos.
  So for now, we check what type the galerkin_state is and from that, we
  set the types for the galerkin residual and jacobian.
  For example, if galerkin_state = eigen_vector, then it would make sense
  to set galerkin_residual = galerkin_state and galerkin_jacobian = eigen_matrix
*/

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

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template<typename T>
struct select_galerkin_types<
  T,
  mpl::enable_if_t<
    ::pressio::is_array_pybind<T>::value
    >
  >
{
  using residual_type = mpl::remove_cvref_t<T>;
  using jacobian_type = mpl::remove_cvref_t<T>;
};
#endif

template<class ...Args> struct ComposeContTimeExplicit;
template<class ...Args> struct ComposeContTimeImplicit;
template<std::size_t, class ...Args> struct ComposeDiscTime;
template<typename ...Args>
using ComposeContTimeExplicit_t = typename ComposeContTimeExplicit<Args...>::type;
template<typename ...Args>
using ComposeContTimeImplicit_t = typename ComposeContTimeImplicit<Args...>::type;
template<std::size_t n, class ...Args>
using ComposeDiscTime_t = typename ComposeDiscTime<n, Args...>::type;

//////////////////////////////////////
/// default
//////////////////////////////////////

/// explicit stepping
template<
  class FomSystemType,
  class DecoderType,
  class GalerkinStateType,
  class FomReferenceState
  >
struct ComposeContTimeExplicit<
  FomSystemType, DecoderType, GalerkinStateType, FomReferenceState
  >
{

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

  using type = ::pressio::rom::galerkin::impl::ProblemExplicit<
    0, FomSystemType, GalerkinStateType, DecoderType>;
};

/// implicit stepping
template<
  class FomSystemType,
  class DecoderType,
  class GalerkinStateType,
  class FomReferenceState
  >
struct ComposeContTimeImplicit<
  FomSystemType, DecoderType, GalerkinStateType, FomReferenceState
  >
{

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

  using type = ::pressio::rom::galerkin::impl::ProblemImplicit<
    1, FomSystemType, GalerkinStateType, galerkin_residual_t, galerkin_jacobian_t, DecoderType>;
};

template<
  std::size_t num_states,
  class FomSystemType,
  class DecoderType,
  class GalerkinStateType,
  class FomReferenceState
  >
struct ComposeDiscTime<
  num_states, FomSystemType, DecoderType, GalerkinStateType, FomReferenceState
  >
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

  using type = ::pressio::rom::galerkin::impl::ProblemImplicit<
    2, FomSystemType, GalerkinStateType,
    galerkin_residual_t, galerkin_jacobian_t, DecoderType,
    ::pressio::ode::StepperTotalNumberOfStates<num_states>>;
};


//////////////////////////////////////
/// masked
//////////////////////////////////////

// explicit stepping
template<
  class FomSystemType,
  class DecoderType,
  class GalerkinStateType,
  class FomReferenceState,
  class MaskerType,
  class ProjectorType
  >
struct ComposeContTimeExplicit<
  FomSystemType, DecoderType, GalerkinStateType, FomReferenceState, MaskerType, ProjectorType
  >
{

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

  using type = ::pressio::rom::galerkin::impl::ProblemExplicit<
    3, FomSystemType, GalerkinStateType, DecoderType, MaskerType, ProjectorType>;
};

/// implicit stepping
template<
  class FomSystemType,
  class DecoderType,
  class GalerkinStateType,
  class FomReferenceState,
  class MaskerType,
  class ProjectorType
  >
struct ComposeContTimeImplicit<
  FomSystemType, DecoderType, GalerkinStateType, FomReferenceState, MaskerType, ProjectorType
  >
{

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

  using type = ::pressio::rom::galerkin::impl::ProblemImplicit<
    4, FomSystemType, GalerkinStateType, galerkin_residual_t,
    galerkin_jacobian_t, DecoderType, MaskerType, ProjectorType>;
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
struct ComposeDiscTime<
  num_states, FomSystemType, DecoderType, GalerkinStateType, FomReferenceState,
  ProjectorType, MaskerType
  >
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

  using type = ::pressio::rom::galerkin::impl::ProblemImplicit<
    5, FomSystemType, GalerkinStateType,
    galerkin_residual_t, galerkin_jacobian_t, DecoderType, ProjectorType, MaskerType,
    ::pressio::ode::StepperTotalNumberOfStates<num_states>>;
};


//////////////////////////////////////
/// hyper-reduced
//////////////////////////////////////

// explicit stepping
template<
  class FomSystemType,
  class DecoderType,
  class GalerkinStateType,
  class FomReferenceState,
  class ProjectorType
  >
struct ComposeContTimeExplicit<
  FomSystemType, DecoderType, GalerkinStateType, FomReferenceState, ProjectorType
  >
{

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

  using type = ::pressio::rom::galerkin::impl::ProblemExplicit<
    6, FomSystemType, GalerkinStateType, DecoderType, ProjectorType>;
};

/// implicit stepping
template<
  class FomSystemType,
  class DecoderType,
  class GalerkinStateType,
  class FomReferenceState,
  class ProjectorType
  >
struct ComposeContTimeImplicit<
    FomSystemType, DecoderType, GalerkinStateType, FomReferenceState, ProjectorType
  >
{

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

  using type = ::pressio::rom::galerkin::impl::ProblemImplicit<
    7, FomSystemType, GalerkinStateType, galerkin_residual_t, galerkin_jacobian_t,
    DecoderType, ProjectorType>;
};

template<
  std::size_t num_states,
  class FomSystemType,
  class DecoderType,
  class GalerkinStateType,
  class FomReferenceState,
  class ProjectorType
  >
struct ComposeDiscTime<
  num_states, FomSystemType, DecoderType, GalerkinStateType,
  FomReferenceState, ProjectorType
  >
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

  using type = ::pressio::rom::galerkin::impl::ProblemImplicit<
    8, FomSystemType, GalerkinStateType,
    galerkin_residual_t, galerkin_jacobian_t, DecoderType, ProjectorType,
    ::pressio::ode::StepperTotalNumberOfStates<num_states>>;
};


void ensure_explicit_or_throw(const std::string & str, ::pressio::ode::StepScheme name)
{
  const auto tmp_b = ::pressio::ode::is_explicit_scheme(name);
  if (!tmp_b){
    throw std::runtime_error(str + " requires an explicit stepper");
  }
}

void ensure_implicit_or_throw(const std::string & str, ::pressio::ode::StepScheme name)
{
  const auto tmp_b = ::pressio::ode::is_explicit_scheme(name);
  if (tmp_b){
    throw std::runtime_error(str + " requires an implicit stepper");
  }
}

}}}}
#endif  // ROM_IMPL_ROM_GALERKIN_COMPOSE_HPP_
