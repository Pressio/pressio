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

#include "../rom_problem_tags.hpp"
#include "../rom_galerkin_types_selector.hpp"

#include "./policies/rom_galerkin_fom_velocity_policy.hpp"
#include "./policies/rom_galerkin_fom_apply_jacobian_policy.hpp"
#include "./policies/rom_galerkin_velocity_policy.hpp"

//traits
#include "./traits/rom_galerkin_common_traits.hpp"
#include "./traits/rom_galerkin_default_problem_explicit_step_traits.hpp"
#include "./traits/rom_galerkin_default_problem_implicit_step_traits.hpp"
#include "./traits/rom_galerkin_hypred_vel_problem_explicit_step_traits.hpp"
#include "./traits/rom_galerkin_hypred_vel_problem_implicit_step_traits.hpp"
#include "./traits/rom_galerkin_masked_vel_problem_explicit_step_traits.hpp"
#include "./traits/rom_galerkin_masked_vel_problem_implicit_step_traits.hpp"
//problems
#include "./rom_galerkin_default_problem_explicit_stepping.hpp"
#include "./rom_galerkin_default_problem_implicit_stepping.hpp"
#include "./rom_galerkin_hypred_vel_problem_explicit_stepping.hpp"
#include "./rom_galerkin_hypred_vel_problem_implicit_stepping.hpp"
#include "./rom_galerkin_masked_vel_problem_explicit_stepping.hpp"
#include "./rom_galerkin_masked_vel_problem_implicit_stepping.hpp"

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{ namespace conttime{

template <typename tag>
struct supported_implicit_stepper_tag{
  static_assert
  (std::is_same<tag, ::pressio::ode::implicitmethods::Euler>::value or
   std::is_same<tag, ::pressio::ode::implicitmethods::BDF2>::value or
   std::is_same<tag, ::pressio::ode::implicitmethods::CrankNicolson>::value,
   "The implicit stepper tag you are passing to create the galerkin problem \
is not supported: this can be because the Galerkin implementation does \
not support it, or because you added a new ode scheme in the ode package \
but forgot to update the list of implicit tags supported by Galerkin which \
currently contains: BDF1, BDF2 or CrankNicolson");
  static constexpr auto value = true;
};

template <typename tag>
struct supported_explicit_stepper_tag{
  static_assert
  (std::is_same<tag, ::pressio::ode::explicitmethods::Euler>::value or
   std::is_same<tag, ::pressio::ode::explicitmethods::RungeKutta4>::value or
   std::is_same<tag, ::pressio::ode::explicitmethods::AdamsBashforth2>::value,
   "The explicit stepper tag you are passing to create the galerkin problem \
is not supported: this can be because the Galerkin implementation does \
not support it, or because you added a new ode scheme in the ode package \
but forgot to update the list of explicit tags supported by Galerkin which \
currently contains: Forward Euler, RK4, AdamsBashforth2");

  static constexpr auto value = true;
};

template <typename tag>
struct valid_stepper_tag{
  static constexpr auto value =
    supported_explicit_stepper_tag<tag>::value or
    supported_implicit_stepper_tag<tag>::value;
};

//------------------------
// COMPOSE
//------------------------
/*
  NOTE (A): pay attention below since for implicit time stepping we need to
  know the type of the reduced residual and jacobian.
  Obviously, this should be something that is compatible
  with the galerkin state and such that we can solve the reduced problem
  easily: typically, the galerkin state is a shared-mem array, e.g. eigen or Kokkos.
  So when the user does not pass them explicitly, we check what type is
  the galerkin_state_t and from that set the types for the galerkin residual and jacobian.
  For example, if galerkin_state = Vector<eigen_vector>, then it would make sense
  to set galerkin_residual = galerkin_state and
  galerkin_jacobian = DenseMatrix<eigen_matrix>
*/

template<
  typename problem_tag,
  typename dummy,
  typename ode_tag,
  typename fom_system_t,
  typename ...Args>
struct compose
{
  /* if we fall here, it means something is wrong
     because it could not match any specialization below.
     Use assertions to tell users what is wrong.
  */
  // 1. check that the ode_tag is a valid stepper tag
  static constexpr auto is_ode_tag =
    ::pressio::ode::predicates::is_stepper_tag<ode_tag>::value;
  static_assert
    (is_ode_tag,
     "Galerkin with continuous-time API: to set the stepping scheme, \
it seems you are using a tag type that is not a valid ode stepper tag.   \
This error is typically caused by the way you create the galerkin problem: \
e.g. createDefaultProblem<ode_tag>(...)");

  using type = void;
};

//*********************************************************
// EXPLICIT TIME STEPPING
//*********************************************************
/****
     Galerkin default, pressio ops, explicit stepping
***/
template<
  typename stepper_tag,
  typename fom_system_type,
  typename decoder_type,
  typename galerkin_state_type
  >
struct compose<
  ::pressio::rom::galerkin::impl::Default,
  mpl::enable_if_t<
    ::pressio::ode::predicates::is_explicit_stepper_tag<stepper_tag>::value
    >,
  stepper_tag, fom_system_type, decoder_type, galerkin_state_type>
{
  static_assert(supported_explicit_stepper_tag<stepper_tag>::value, "");

  using type =
    ::pressio::rom::galerkin::impl::DefaultProblemExplicitStepContinuousTimeApi<
    stepper_tag, fom_system_type, galerkin_state_type, decoder_type, void>;
};

/****
     Galerkin default, user-defined ops, explicit stepping
***/
template<
  typename stepper_tag,
  typename fom_system_type,
  typename decoder_type,
  typename galerkin_state_type,
  typename ud_ops_type
  >
struct compose<
  ::pressio::rom::galerkin::impl::Default,
  mpl::enable_if_t<
    ::pressio::ode::predicates::is_explicit_stepper_tag<stepper_tag>::value
    >,
  stepper_tag, fom_system_type, decoder_type, galerkin_state_type, ud_ops_type>
{
  static_assert(supported_explicit_stepper_tag<stepper_tag>::value, "");
  static_assert(mpl::not_void<ud_ops_type>::value, "ud_ops_type cannot be void");

  using type =
    ::pressio::rom::galerkin::impl::DefaultProblemExplicitStepContinuousTimeApi<
    stepper_tag, fom_system_type, galerkin_state_type, decoder_type, ud_ops_type>;
};

/****
     hyperReducedVelocity Galerkin, pressio ops, explicit stepping
***/
template<
  typename stepper_tag,
  typename fom_system_type,
  typename decoder_type,
  typename galerkin_state_type,
  typename projector_type
  >
struct compose<
  ::pressio::rom::galerkin::impl::HyperReducedVelocity,
  mpl::enable_if_t<
    ::pressio::ode::predicates::is_explicit_stepper_tag<stepper_tag>::value
    >,
  stepper_tag, fom_system_type, decoder_type, galerkin_state_type, projector_type>
{
  static_assert(supported_explicit_stepper_tag<stepper_tag>::value, "");

  using type =
    ::pressio::rom::galerkin::impl::HypRedVeloProblemExplicitStepContinuousTimeApi<
    stepper_tag, fom_system_type, galerkin_state_type, decoder_type, projector_type, void>;
};

/****
     maskedVelocity Galerkin, pressio ops, explicit stepping
***/
template<
  typename stepper_tag,
  typename fom_system_type,
  typename decoder_type,
  typename galerkin_state_type,
  typename masker_type,
  typename projector_type
  >
struct compose<
  ::pressio::rom::galerkin::impl::MaskedVelocity,
  mpl::enable_if_t<
    ::pressio::ode::predicates::is_explicit_stepper_tag<stepper_tag>::value
    >,
  stepper_tag, fom_system_type, decoder_type, galerkin_state_type, masker_type, projector_type>
{
  static_assert(supported_explicit_stepper_tag<stepper_tag>::value, "");

  using type =
    ::pressio::rom::galerkin::impl::MaskedVeloProblemExplicitStepContinuousTimeApi<
    stepper_tag, fom_system_type, galerkin_state_type, decoder_type,
    masker_type, projector_type, void>;
};

//*********************************************************
// IMPLICIT TIME STEPPING
//*********************************************************
/****
     Galerkin default, pressio ops, implicit stepping
***/
template<
  typename stepper_tag,
  typename fom_system_type,
  typename decoder_type,
  typename galerkin_state_type
  >
struct compose<
  ::pressio::rom::galerkin::impl::Default,
  mpl::enable_if_t<
    ::pressio::ode::predicates::is_implicit_stepper_tag<stepper_tag>::value
    >,
  stepper_tag, fom_system_type, decoder_type, galerkin_state_type>
{
  static_assert(supported_implicit_stepper_tag<stepper_tag>::value, "");

  // infer residual and jacobian from state, see NOTE (A) at top
  using galerkin_residual_t =
    typename ::pressio::rom::galerkin::impl::select_galerkin_types<galerkin_state_type>::residual_type;
  using galerkin_jacobian_t =
    typename ::pressio::rom::galerkin::impl::select_galerkin_types<galerkin_state_type>::jacobian_type;

  using type =
    ::pressio::rom::galerkin::impl::DefaultProblemImplicitStepContinuousTimeApi<
    stepper_tag, fom_system_type,
    galerkin_state_type, galerkin_residual_t, galerkin_jacobian_t, decoder_type, void>;
};

/****
     Galerkin default, user-defined ops, implicit stepping
***/
template<
  typename stepper_tag,
  typename fom_system_type,
  typename decoder_type,
  typename galerkin_state_type,
  typename ud_ops_type
  >
struct compose<
  ::pressio::rom::galerkin::impl::Default,
  mpl::enable_if_t<
    ::pressio::ode::predicates::is_implicit_stepper_tag<stepper_tag>::value
    >,
  stepper_tag, fom_system_type, decoder_type, galerkin_state_type, ud_ops_type>
{
  static_assert(supported_implicit_stepper_tag<stepper_tag>::value, "");
  static_assert(mpl::not_void<ud_ops_type>::value, "ud_ops_type cannot be void");

  // infer residual and jacobian from state, see NOTE (A) at top
  using galerkin_residual_t =
    typename ::pressio::rom::galerkin::impl::select_galerkin_types<galerkin_state_type>::residual_type;
  using galerkin_jacobian_t =
    typename ::pressio::rom::galerkin::impl::select_galerkin_types<galerkin_state_type>::jacobian_type;

  using type =
    ::pressio::rom::galerkin::impl::DefaultProblemImplicitStepContinuousTimeApi<
    stepper_tag, fom_system_type,
    galerkin_state_type, galerkin_residual_t, galerkin_jacobian_t,
    decoder_type, ud_ops_type>;
};

/****
     hyperReducedVelocity Galerkin, pressio ops, implicit stepping
***/
template<
  typename stepper_tag,
  typename fom_system_type,
  typename decoder_type,
  typename galerkin_state_type,
  typename projector_type
  >
struct compose<
  ::pressio::rom::galerkin::impl::HyperReducedVelocity,
  mpl::enable_if_t<
    ::pressio::ode::predicates::is_implicit_stepper_tag<stepper_tag>::value
    >,
  stepper_tag, fom_system_type, decoder_type, galerkin_state_type, projector_type>
{
  static_assert(supported_implicit_stepper_tag<stepper_tag>::value, "Invalid stepper tag");

  using galerkin_residual_t =
    typename ::pressio::rom::galerkin::impl::select_galerkin_types<galerkin_state_type>::residual_type;
  using galerkin_jacobian_t =
    typename ::pressio::rom::galerkin::impl::select_galerkin_types<galerkin_state_type>::jacobian_type;

  using type =
    ::pressio::rom::galerkin::impl::HypRedVeloProblemImplicitStepContinuousTimeApi<
    stepper_tag, fom_system_type, galerkin_state_type, galerkin_residual_t,
    galerkin_jacobian_t, decoder_type, projector_type, void>;
};

/****
     hyperReducedVelocity Galerkin, user-defined ops, implicit stepping
***/
template<
  typename stepper_tag,
  typename fom_system_type,
  typename decoder_type,
  typename galerkin_state_type,
  typename projector_type,
  typename ud_ops_type
  >
struct compose<
  ::pressio::rom::galerkin::impl::HyperReducedVelocity,
  mpl::enable_if_t<
    ::pressio::ode::predicates::is_implicit_stepper_tag<stepper_tag>::value
    >,
  stepper_tag, fom_system_type, decoder_type, galerkin_state_type, projector_type, ud_ops_type>
{
  static_assert(supported_implicit_stepper_tag<stepper_tag>::value, "Invalid stepper tag");
  static_assert(mpl::not_void<ud_ops_type>::value, "ud_ops_type cannot be void");

  // infer residual and jacobian from state, see NOTE (A) at top
  using galerkin_residual_t =
    typename ::pressio::rom::galerkin::impl::select_galerkin_types<galerkin_state_type>::residual_type;
  using galerkin_jacobian_t =
    typename ::pressio::rom::galerkin::impl::select_galerkin_types<galerkin_state_type>::jacobian_type;

  using type =
    ::pressio::rom::galerkin::impl::HypRedVeloProblemImplicitStepContinuousTimeApi<
    stepper_tag, fom_system_type, galerkin_state_type, galerkin_residual_t,
    galerkin_jacobian_t, decoder_type, projector_type, ud_ops_type>;
};

/****
     maskedVelocity Galerkin, pressio ops, implicit stepping
***/
template<
  typename stepper_tag,
  typename fom_system_type,
  typename decoder_type,
  typename galerkin_state_type,
  typename masker_type,
  typename projector_type
  >
struct compose<
  ::pressio::rom::galerkin::impl::MaskedVelocity,
  mpl::enable_if_t<
    ::pressio::ode::predicates::is_implicit_stepper_tag<stepper_tag>::value
    >,
  stepper_tag, fom_system_type, decoder_type, galerkin_state_type, masker_type, projector_type>
{
  static_assert(supported_implicit_stepper_tag<stepper_tag>::value, "Invalid stepper tag");

  // infer residual and jacobian from state, see NOTE (A) at top
  using galerkin_residual_t =
    typename ::pressio::rom::galerkin::impl::select_galerkin_types<galerkin_state_type>::residual_type;
  using galerkin_jacobian_t =
    typename ::pressio::rom::galerkin::impl::select_galerkin_types<galerkin_state_type>::jacobian_type;

  using type =
    ::pressio::rom::galerkin::impl::MaskedVeloProblemImplicitStepContinuousTimeApi<
    stepper_tag, fom_system_type, galerkin_state_type, galerkin_residual_t,
    galerkin_jacobian_t, decoder_type, masker_type, projector_type, void>;
};

//-------------------------------------------------------
} // end namespace pressio::rom::galerkin::impl::conttime
//-------------------------------------------------------

// default continuous-time
template<typename ...Args>
using composeDefaultProblemContTime =
  impl::conttime::compose<
  impl::Default, void,
  typename std::remove_cv<typename std::remove_reference<Args>::type>::type...>;

template<typename ...Args>
using composeDefaultProblemContTime_t =
  typename composeDefaultProblemContTime<Args...>::type;


// hr velocity continuous-time
template<typename ...Args>
using composeHyperReducedVelocityProblemContTime =
  impl::conttime::compose<
  impl::HyperReducedVelocity, void,
  typename std::remove_cv<typename std::remove_reference<Args>::type>::type...>;

template<typename ...Args>
using composeHyperReducedVelocityProblemContTime_t =
  typename composeHyperReducedVelocityProblemContTime<Args...>::type;


// masked velocity continuous-time
template<typename ...Args>
using composeMaskedVelocityProblemContTime =
  impl::conttime::compose<
  impl::MaskedVelocity, void,
  typename std::remove_cv<typename std::remove_reference<Args>::type>::type...>;

template<typename ...Args>
using composeMaskedVelocityProblemContTime_t =
  typename composeMaskedVelocityProblemContTime<Args...>::type;

}}}} // end namespace pressio::rom::galerkin::impl
#endif  // ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_ROM_COMPOSE_IMPL_HPP_
