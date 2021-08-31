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

#ifndef ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_ROM_COMPOSE_IMPL_HPP_
#define ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_ROM_COMPOSE_IMPL_HPP_

#include "./rom_problem_members_common_mixins.hpp"
#include "./rom_galerkin_problem_members.hpp"

#include "./rom_galerkin_decorator_projected.hpp"
#include "./rom_galerkin_decorator_masked.hpp"
#include "./rom_galerkin_fom_velocity_policy.hpp"
#include "./rom_galerkin_fom_apply_jacobian_policy.hpp"
#include "./rom_galerkin_velocity_policy.hpp"
#include "./rom_galerkin_policy_residual.hpp"
#include "./rom_galerkin_policy_jacobian.hpp"

#include "./rom_galerkin_cont_time_problem_traits.hpp"
#include "./rom_galerkin_cont_time_problem.hpp"

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

//////////////////////////////////////
/// Cont-time, Galerkin default
//////////////////////////////////////
template<class stepper_tag, class ...Args>
struct ComposeDefaultProblemContTime;

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

  using type = ::pressio::rom::galerkin::impl::ProblemContinuousTimeApi<
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

  using type = ::pressio::rom::galerkin::impl::ProblemContinuousTimeApi<
    1, StepperTag, FomSystemType, GalerkinStateType,
    galerkin_residual_t, galerkin_jacobian_t, DecoderType>;
};

//////////////////////////////////////
/// Cont-time, masked Galerkin
//////////////////////////////////////
template<class stepper_tag, class ...Args>
struct ComposeMaskedProblemContTime;

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

  using type = ::pressio::rom::galerkin::impl::ProblemContinuousTimeApi<
    2, StepperTag, FomSystemType, GalerkinStateType, DecoderType, MaskerType, ProjectorType>;
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


  using type = ::pressio::rom::galerkin::impl::ProblemContinuousTimeApi<
    3, StepperTag, FomSystemType, GalerkinStateType,
    galerkin_residual_t, galerkin_jacobian_t, DecoderType,
    MaskerType, ProjectorType>;
};


//////////////////////////////////////
/// Cont-time, hyper-reduced Galerkin
//////////////////////////////////////
template<class stepper_tag, class ...Args>
struct ComposeHypRedVeloProblemContTime;

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

  using type = ::pressio::rom::galerkin::impl::ProblemContinuousTimeApi<
    4, StepperTag, FomSystemType, GalerkinStateType, DecoderType, ProjectorType>;
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

  using type = ::pressio::rom::galerkin::impl::ProblemContinuousTimeApi<
    5, StepperTag, FomSystemType, GalerkinStateType,
    galerkin_residual_t, galerkin_jacobian_t, DecoderType, ProjectorType>;
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

}}}
#endif  // ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_ROM_COMPOSE_IMPL_HPP_




//namespace conttime{

// template <typename tag>
// struct valid_stepper_tag{
//   static constexpr auto value =
//     supported_explicit_stepper_tag<tag>::value or
//     supported_implicit_stepper_tag<tag>::value;
// };

// //------------------------
// // COMPOSE
// //------------------------
// /*
//   NOTE (A): pay attention below since for implicit time stepping we need to
//   know the type of the reduced residual and jacobian.
//   Obviously, this should be something that is compatible
//   with the galerkin state and such that we can solve the reduced problem
//   easily: typically, the galerkin state is a shared-mem array, e.g. eigen or Kokkos.
//   So when the user does not pass them explicitly, we check what type is
//   the galerkin_state_t and from that set the types for the galerkin residual and jacobian.
//   For example, if galerkin_state = Vector<eigen_vector>, then it would make sense
//   to set galerkin_residual = galerkin_state and
//   galerkin_jacobian = DenseMatrix<eigen_matrix>
// */

// template<
//   typename problem_tag,
//   typename dummy,
//   typename ode_tag,
//   typename fom_system_t,
//   typename ...Args>
// struct compose
// {
//   /* if we fall here, it means something is wrong
//      because it could not match any specialization below.
//      Use assertions to tell users what is wrong.
//   */
//   // 1. check that the ode_tag is a valid stepper tag
//   static constexpr auto is_ode_tag =
//     ::pressio::ode::predicates::is_stepper_tag<ode_tag>::value;
//   static_assert
//     (is_ode_tag,
//      "Galerkin with continuous-time API: to set the stepping scheme, it seems you are using a tag type that is not a valid ode stepper tag. This error is typically caused by the way you create the galerkin problem: e.g. createDefaultProblem<ode_tag>(...)");
//   using type = void;
// };


// //*********************************************************
// // EXPLICIT TIME STEPPING
// //*********************************************************
// /****
//      Galerkin default, explicit stepping
// ***/
// template<
//   typename stepper_tag,
//   typename fom_system_type,
//   typename decoder_type,
//   typename galerkin_state_type
//   >
// struct compose<
//   ::pressio::rom::galerkin::impl::Default,
//   mpl::enable_if_t<
//     ::pressio::ode::is_explicit_stepper_tag<stepper_tag>::value
//     >,
//   stepper_tag, fom_system_type, decoder_type, galerkin_state_type>
// {
//   static_assert(supported_explicit_stepper_tag<stepper_tag>::value, "");

//   using type =
//     ::pressio::rom::galerkin::impl::DefaultProblemExplicitStepContinuousTimeApi<
//     stepper_tag, fom_system_type, galerkin_state_type, decoder_type, void>;
// };

// /****
//      hyperReducedVelocity Galerkin, pressio ops, explicit stepping
// ***/
// template<
//   typename stepper_tag,
//   typename fom_system_type,
//   typename decoder_type,
//   typename galerkin_state_type,
//   typename projector_type
//   >
// struct compose<
//   ::pressio::rom::galerkin::impl::HyperReducedVelocity,
//   mpl::enable_if_t<
//     ::pressio::ode::is_explicit_stepper_tag<stepper_tag>::value
//     >,
//   stepper_tag, fom_system_type, decoder_type, galerkin_state_type, projector_type>
// {
//   static_assert(supported_explicit_stepper_tag<stepper_tag>::value, "");

//   using type =
//     ::pressio::rom::galerkin::impl::HypRedVeloProblemExplicitStepContinuousTimeApi<
//     stepper_tag, fom_system_type, galerkin_state_type, decoder_type, projector_type, void>;
// };

// /****
//      maskedVelocity Galerkin, pressio ops, explicit stepping
// ***/
// template<
//   typename stepper_tag,
//   typename fom_system_type,
//   typename decoder_type,
//   typename galerkin_state_type,
//   typename masker_type,
//   typename projector_type
//   >
// struct compose<
//   ::pressio::rom::galerkin::impl::MaskedVelocity,
//   mpl::enable_if_t<
//     ::pressio::ode::is_explicit_stepper_tag<stepper_tag>::value
//     >,
//   stepper_tag, fom_system_type, decoder_type, galerkin_state_type, masker_type, projector_type>
// {
//   static_assert(supported_explicit_stepper_tag<stepper_tag>::value, "");

//   using type =
//     ::pressio::rom::galerkin::impl::MaskedVeloProblemExplicitStepContinuousTimeApi<
//     stepper_tag, fom_system_type, galerkin_state_type, decoder_type,
//     masker_type, projector_type, void>;
// };

// //*********************************************************
// // IMPLICIT TIME STEPPING
// //*********************************************************
// /****
//      Galerkin default, pressio ops, implicit stepping
// ***/
// template<
//   typename stepper_tag,
//   typename fom_system_type,
//   typename decoder_type,
//   typename galerkin_state_type
//   >
// struct compose<
//   ::pressio::rom::galerkin::impl::Default,
//   mpl::enable_if_t<
//     ::pressio::ode::predicates::is_implicit_stepper_tag<stepper_tag>::value
//     >,
//   stepper_tag, fom_system_type, decoder_type, galerkin_state_type>
// {
//   static_assert(supported_implicit_stepper_tag<stepper_tag>::value, "");

//   // infer residual and jacobian from state, see NOTE (A) at top
//   using galerkin_residual_t =
//     typename ::pressio::rom::galerkin::impl::select_galerkin_types<galerkin_state_type>::residual_type;
//   using galerkin_jacobian_t =
//     typename ::pressio::rom::galerkin::impl::select_galerkin_types<galerkin_state_type>::jacobian_type;

//   using type =
//     ::pressio::rom::galerkin::impl::DefaultProblemImplicitStepContinuousTimeApi<
//     stepper_tag, fom_system_type,
//     galerkin_state_type, galerkin_residual_t, galerkin_jacobian_t, decoder_type, void>;
// };

// /****
//      Galerkin default, user-defined ops, implicit stepping
// ***/
// template<
//   typename stepper_tag,
//   typename fom_system_type,
//   typename decoder_type,
//   typename galerkin_state_type,
//   typename ud_ops_type
//   >
// struct compose<
//   ::pressio::rom::galerkin::impl::Default,
//   mpl::enable_if_t<
//     ::pressio::ode::predicates::is_implicit_stepper_tag<stepper_tag>::value
//     >,
//   stepper_tag, fom_system_type, decoder_type, galerkin_state_type, ud_ops_type>
// {
//   static_assert(supported_implicit_stepper_tag<stepper_tag>::value, "");
//   static_assert(mpl::not_void<ud_ops_type>::value, "ud_ops_type cannot be void");

//   // infer residual and jacobian from state, see NOTE (A) at top
//   using galerkin_residual_t =
//     typename ::pressio::rom::galerkin::impl::select_galerkin_types<galerkin_state_type>::residual_type;
//   using galerkin_jacobian_t =
//     typename ::pressio::rom::galerkin::impl::select_galerkin_types<galerkin_state_type>::jacobian_type;

//   using type =
//     ::pressio::rom::galerkin::impl::DefaultProblemImplicitStepContinuousTimeApi<
//     stepper_tag, fom_system_type,
//     galerkin_state_type, galerkin_residual_t, galerkin_jacobian_t,
//     decoder_type, ud_ops_type>;
// };

// /****
//      hyperReducedVelocity Galerkin, pressio ops, implicit stepping
// ***/
// template<
//   typename stepper_tag,
//   typename fom_system_type,
//   typename decoder_type,
//   typename galerkin_state_type,
//   typename projector_type
//   >
// struct compose<
//   ::pressio::rom::galerkin::impl::HyperReducedVelocity,
//   mpl::enable_if_t<
//     ::pressio::ode::predicates::is_implicit_stepper_tag<stepper_tag>::value
//     >,
//   stepper_tag, fom_system_type, decoder_type, galerkin_state_type, projector_type>
// {
//   static_assert(supported_implicit_stepper_tag<stepper_tag>::value, "Invalid stepper tag");

//   using galerkin_residual_t =
//     typename ::pressio::rom::galerkin::impl::select_galerkin_types<galerkin_state_type>::residual_type;
//   using galerkin_jacobian_t =
//     typename ::pressio::rom::galerkin::impl::select_galerkin_types<galerkin_state_type>::jacobian_type;

//   using type =
//     ::pressio::rom::galerkin::impl::HypRedVeloProblemImplicitStepContinuousTimeApi<
//     stepper_tag, fom_system_type, galerkin_state_type, galerkin_residual_t,
//     galerkin_jacobian_t, decoder_type, projector_type, void>;
// };

// /****
//      hyperReducedVelocity Galerkin, user-defined ops, implicit stepping
// ***/
// template<
//   typename stepper_tag,
//   typename fom_system_type,
//   typename decoder_type,
//   typename galerkin_state_type,
//   typename projector_type,
//   typename ud_ops_type
//   >
// struct compose<
//   ::pressio::rom::galerkin::impl::HyperReducedVelocity,
//   mpl::enable_if_t<
//     ::pressio::ode::predicates::is_implicit_stepper_tag<stepper_tag>::value
//     >,
//   stepper_tag, fom_system_type, decoder_type, galerkin_state_type, projector_type, ud_ops_type>
// {
//   static_assert(supported_implicit_stepper_tag<stepper_tag>::value, "Invalid stepper tag");
//   static_assert(mpl::not_void<ud_ops_type>::value, "ud_ops_type cannot be void");

//   // infer residual and jacobian from state, see NOTE (A) at top
//   using galerkin_residual_t =
//     typename ::pressio::rom::galerkin::impl::select_galerkin_types<galerkin_state_type>::residual_type;
//   using galerkin_jacobian_t =
//     typename ::pressio::rom::galerkin::impl::select_galerkin_types<galerkin_state_type>::jacobian_type;

//   using type =
//     ::pressio::rom::galerkin::impl::HypRedVeloProblemImplicitStepContinuousTimeApi<
//     stepper_tag, fom_system_type, galerkin_state_type, galerkin_residual_t,
//     galerkin_jacobian_t, decoder_type, projector_type, ud_ops_type>;
// };

// /****
//      maskedVelocity Galerkin, pressio ops, implicit stepping
// ***/
// template<
//   typename stepper_tag,
//   typename fom_system_type,
//   typename decoder_type,
//   typename galerkin_state_type,
//   typename masker_type,
//   typename projector_type
//   >
// struct compose<
//   ::pressio::rom::galerkin::impl::MaskedVelocity,
//   mpl::enable_if_t<
//     ::pressio::ode::predicates::is_implicit_stepper_tag<stepper_tag>::value
//     >,
//   stepper_tag, fom_system_type, decoder_type, galerkin_state_type, masker_type, projector_type>
// {
//   static_assert(supported_implicit_stepper_tag<stepper_tag>::value, "Invalid stepper tag");

//   // infer residual and jacobian from state, see NOTE (A) at top
//   using galerkin_residual_t =
//     typename ::pressio::rom::galerkin::impl::select_galerkin_types<galerkin_state_type>::residual_type;
//   using galerkin_jacobian_t =
//     typename ::pressio::rom::galerkin::impl::select_galerkin_types<galerkin_state_type>::jacobian_type;

//   using type =
//     ::pressio::rom::galerkin::impl::MaskedVeloProblemImplicitStepContinuousTimeApi<
//     stepper_tag, fom_system_type, galerkin_state_type, galerkin_residual_t,
//     galerkin_jacobian_t, decoder_type, masker_type, projector_type, void>;
// };

//-------------------------------------------------------
} // end namespace pressio::rom::galerkin::impl::conttime
//-------------------------------------------------------

// // default continuous-time
// template<typename ...Args>
// using composeDefaultProblemContTime =
//   impl::conttime::compose<
//     impl::Default, void,
//     typename std::remove_cv<typename std::remove_reference<Args>::type>::type...
//   >;

// template<typename ...Args>
// using composeDefaultProblemContTime_t =
//   typename composeDefaultProblemContTime<Args...>::type;

// // hr velocity continuous-time
// template<typename ...Args>
// using composeHyperReducedVelocityProblemContTime =
//   impl::conttime::compose<
//   impl::HyperReducedVelocity, void,
//   typename std::remove_cv<typename std::remove_reference<Args>::type>::type...>;

// template<typename ...Args>
// using composeHyperReducedVelocityProblemContTime_t =
//   typename composeHyperReducedVelocityProblemContTime<Args...>::type;

// // masked velocity continuous-time
// template<typename ...Args>
// using composeMaskedVelocityProblemContTime =
//   impl::conttime::compose<
//   impl::MaskedVelocity, void,
//   typename std::remove_cv<typename std::remove_reference<Args>::type>::type...>;

// template<typename ...Args>
// using composeMaskedVelocityProblemContTime_t =
//   typename composeMaskedVelocityProblemContTime<Args...>::type;
//}// end namespace conttime

// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
// template<typename T>
// struct select_galerkin_types<
//   T,
//   mpl::enable_if_t<
//     ::pressio::containers::predicates::is_rank1_tensor_wrapper_pybind<T>::value
//     >
//   >
// {
//   using scalar_t = typename ::pressio::containers::details::traits<T>::scalar_t;
//   // for now use residual_type = state_type
//   using residual_type = T;
//   // the galerkin jacobian is a pybind11 tensor column-major
//   using native_j_t = pybind11::array_t<scalar_t, pybind11::array::f_style>;
//   using jacobian_type = ::pressio::containers::Tensor<2, native_j_t>;
// };
// #endif
