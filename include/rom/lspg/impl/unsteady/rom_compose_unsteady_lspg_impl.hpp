/*
//@HEADER
// ************************************************************************
//
// rom_compose_unsteady_lspg_impl.hpp
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

#ifndef ROM_LSPG_IMPL_UNSTEADY_ROM_COMPOSE_UNSTEADY_LSPG_IMPL_HPP_
#define ROM_LSPG_IMPL_UNSTEADY_ROM_COMPOSE_UNSTEADY_LSPG_IMPL_HPP_

#include "../rom_problem_tags.hpp"

//----------------------------------------
// include for continuous_time_api
//----------------------------------------
#include "./continuous_time_api/discrete_time_functions/rom_lspg_time_discrete_residual.hpp"
#include "./continuous_time_api/discrete_time_functions/rom_lspg_time_discrete_jacobian.hpp"
#include "./continuous_time_api/policies/rom_lspg_unsteady_residual_policy_continuous_time_api.hpp"
#include "./continuous_time_api/policies/rom_lspg_unsteady_jacobian_policy_continuous_time_api.hpp"
#include "./continuous_time_api/policies/rom_lspg_unsteady_hyper_reduced_residual_policy_continuous_time_api.hpp"
#include "./continuous_time_api/policies/rom_lspg_unsteady_hyper_reduced_jacobian_policy_continuous_time_api.hpp"
//traits
#include "./continuous_time_api/traits/rom_lspg_unsteady_common_traits_continuous_time_api.hpp"
#include "./continuous_time_api/traits/rom_lspg_unsteady_default_problem_traits_continuous_time_api.hpp"
#include "./continuous_time_api/traits/rom_lspg_unsteady_preconditioned_problem_traits_continuous_time_api.hpp"
#include "./continuous_time_api/traits/rom_lspg_unsteady_masked_problem_traits_continuous_time_api.hpp"
#include "./continuous_time_api/traits/rom_lspg_unsteady_hyper_reduced_problem_traits_continuous_time_api.hpp"
#include "./continuous_time_api/traits/rom_lspg_unsteady_preconditioned_hyper_reduced_problem_traits_continuous_time_api.hpp"
//problems
#include "./continuous_time_api/rom_lspg_unsteady_default_problem_continuous_time_api.hpp"
#include "./continuous_time_api/rom_lspg_unsteady_preconditioned_problem_continuous_time_api.hpp"
#include "./continuous_time_api/rom_lspg_unsteady_masked_problem_continuous_time_api.hpp"
#include "./continuous_time_api/rom_lspg_unsteady_hyper_reduced_problem_continuous_time_api.hpp"
#include "./continuous_time_api/rom_lspg_unsteady_preconditioned_hyper_reduced_problem_continuous_time_api.hpp"

//----------------------------------------
// include for discrete_time_api
//----------------------------------------
#include "./discrete_time_api/policies/rom_lspg_unsteady_residual_policy_discrete_time_api.hpp"
#include "./discrete_time_api/policies/rom_lspg_unsteady_jacobian_policy_discrete_time_api.hpp"
//traits
#include "./discrete_time_api/traits/rom_lspg_unsteady_common_traits_discrete_time_api.hpp"
#include "./discrete_time_api/traits/rom_lspg_unsteady_default_problem_traits_discrete_time_api.hpp"
#include "./discrete_time_api/traits/rom_lspg_unsteady_preconditioned_problem_traits_discrete_time_api.hpp"
#include "./discrete_time_api/traits/rom_lspg_unsteady_masked_problem_traits_discrete_time_api.hpp"
//problems
#include "./discrete_time_api/rom_lspg_unsteady_default_problem_discrete_time_api.hpp"
#include "./discrete_time_api/rom_lspg_unsteady_preconditioned_problem_discrete_time_api.hpp"
#include "./discrete_time_api/rom_lspg_unsteady_masked_problem_discrete_time_api.hpp"


namespace pressio{ namespace rom{ namespace lspg{ namespace impl{

template <typename tag>
struct valid_stepper_tag_continuous_time_api
{
  static_assert
  (std::is_same<tag, ::pressio::ode::implicitmethods::Euler>::value or
   std::is_same<tag, ::pressio::ode::implicitmethods::BDF2>::value or
   std::is_same<tag, ::pressio::ode::implicitmethods::CrankNicolson>::value,
   "The implicit stepper tag you are passing to create the LSPG problem \
is not supported: this can be because the current LSPG implementation does \
not support it, or because you added a new ode scheme in the ode package \
but forgot to update the list of implicit tags supported by LSPG which \
currently contains: BDF1, BDF2 or CrankNicolson");

  static constexpr auto value = true;
};

//##############
//##############

template<
  bool is_probably_discrete_time,
  typename fom_system_t,
  typename apply_jac_operand_t
  >
struct FindAdapterMistakesUnsteady;

template<typename fom_system_t, typename apply_jac_operand_t>
struct FindAdapterMistakesUnsteady<true, fom_system_t, apply_jac_operand_t>
{
  static constexpr bool value =
    (::pressio::rom::why_not_discrete_time_system_with_user_provided_apply_jacobian
     <fom_system_t, apply_jac_operand_t>::value, "");
};

template<typename fom_system_t, typename apply_jac_operand_t>
struct FindAdapterMistakesUnsteady<false, fom_system_t, apply_jac_operand_t>
{
  static constexpr bool value =
    (::pressio::rom::why_not_continuous_time_system_with_user_provided_apply_jacobian
     <fom_system_t, apply_jac_operand_t>::value, "");
};

template<
  typename problem_tag,
  typename dummy,
  typename ode_tag,
  typename fom_system_t,
  typename decoder_type,
  typename ...Args>
struct composeUnsteady
{
  //if we are here, something is wrong, find out what it is

  // check if the stepper tag is wrong
  static_assert
  (::pressio::ode::predicates::is_stepper_tag<ode_tag>::value,
   "\nThe unsteady LSPG problem you are trying to create cannot be created because \
the first template argument is not a valid stepper tag.");

  /*if here, it means that the adapter class is the problem.
    There are two scenarios: either the user wants to use the discrete-time
    or the continuous-time API. But we have no way to know here which one
    they WANT to use and want to print error messages detailed enough.
    So to do that we guess what API they want
    by discriminating on whether the adapter has a discrete_time_residual_type.
    We choose this because it is simple enough so we hope the user
    does not get it wrong!
    We believe that if the system class has
    the discrete_time_residual_type typedef,
    then it is very likely the user is trying to use a discrete-time API.
    If not, then they are trying to use the continuous-time API.
    If the user if actually using the discrete-time API but
    mispells or forgets the "discrete_time_residual_type" typedef,
    then the logic here breaks, and we might need to add logic to handle the
    case if that typedef is not found, because it could mean the user is
    using the continuous-time api or mispelled the typedef or forgot it.
  */
  using error =
    typename std::conditional<
    ::pressio::ode::predicates::has_discrete_time_residual_typedef<fom_system_t>::value,
    FindAdapterMistakesUnsteady<true,  fom_system_t, typename decoder_type::jacobian_type>,
    FindAdapterMistakesUnsteady<false, fom_system_t, typename decoder_type::jacobian_type>
    >::type;

  // I need this static assert otherwise it won't trigger errors
  static_assert(error::value,"");

  // if we get here set to void, but we should never get here in thory
  // because the asserts above should be triggered before
  using type = void;
};

///////////////////////////////////////////////////////
//////       CONTINUOUS TIME API
///////////////////////////////////////////////////////

/***
  default lspg pressio ops
***/
template<
  typename stepper_tag,
  typename fom_system_type,
  typename decoder_type,
  typename lspg_state_type
  >
struct composeUnsteady<
  ::pressio::rom::lspg::impl::Default,
  mpl::enable_if_t<
    ::pressio::rom::constraints::continuous_time_system_with_user_provided_apply_jacobian<
      fom_system_type, typename decoder_type::jacobian_type>::value
    >,
  stepper_tag, fom_system_type, decoder_type, lspg_state_type>
{
  static_assert(valid_stepper_tag_continuous_time_api<stepper_tag>::value,"");

  using type =
    ::pressio::rom::lspg::impl::unsteady::DefaultProblemContinuousTimeApi<
    stepper_tag, fom_system_type, lspg_state_type, decoder_type, void>;
};

/***
  default lspg with user-defined ops
***/
template<
  typename stepper_tag,
  typename fom_system_type,
  typename decoder_type,
  typename lspg_state_type,
  typename ud_ops_type
  >
struct composeUnsteady<
  ::pressio::rom::lspg::impl::Default,
  mpl::enable_if_t<
    ::pressio::rom::constraints::continuous_time_system_with_user_provided_apply_jacobian<
      fom_system_type, typename decoder_type::jacobian_type>::value
    >,
  stepper_tag, fom_system_type, decoder_type, lspg_state_type, ud_ops_type>
{
  static_assert(valid_stepper_tag_continuous_time_api<stepper_tag>::value,"");

  using type =
    ::pressio::rom::lspg::impl::unsteady::DefaultProblemContinuousTimeApi<
    stepper_tag, fom_system_type, lspg_state_type, decoder_type, ud_ops_type>;
};

/***
    preconditioned default lspg pressio ops
***/
template<
  typename stepper_tag,
  typename fom_system_type,
  typename decoder_type,
  typename lspg_state_type,
  typename precond_type
  >
struct composeUnsteady<
  ::pressio::rom::lspg::impl::Preconditioned,
  mpl::enable_if_t<
    ::pressio::rom::constraints::continuous_time_system_with_user_provided_apply_jacobian<
      fom_system_type, typename decoder_type::jacobian_type>::value
    >,
  stepper_tag, fom_system_type, decoder_type, lspg_state_type, precond_type>
{
  static_assert(valid_stepper_tag_continuous_time_api<stepper_tag>::value,"");

  using type =
    ::pressio::rom::lspg::impl::unsteady::PreconditionedProblemContinuousTimeApi<
    stepper_tag, fom_system_type, lspg_state_type, decoder_type, precond_type, void>;
};

/***
    masked lspg pressio ops
***/
template<
  typename stepper_tag,
  typename fom_system_type,
  typename decoder_type,
  typename lspg_state_type,
  typename masker_type
  >
struct composeUnsteady<
  ::pressio::rom::lspg::impl::Masked,
  mpl::enable_if_t<
    ::pressio::rom::constraints::continuous_time_system_with_user_provided_apply_jacobian<
      fom_system_type, typename decoder_type::jacobian_type>::value
    >,
  stepper_tag, fom_system_type, decoder_type, lspg_state_type, masker_type>
{
  static_assert(valid_stepper_tag_continuous_time_api<stepper_tag>::value,"");

  using type =
    ::pressio::rom::lspg::impl::unsteady::MaskedProblemContinuousTimeApi<
    stepper_tag, fom_system_type, lspg_state_type, decoder_type, masker_type, void>;
};


#if defined PRESSIO_ENABLE_TPL_EIGEN or defined PRESSIO_ENABLE_TPL_PYBIND11
/***
    hyper-reduced lspg with pressio ops and provided mapping from sample to stencil mesh.
    currently only enabled for shared-mem data structures: eigen and pybind11
***/
template<
  typename stepper_tag,
  typename fom_system_type,
  typename decoder_type,
  typename lspg_state_type,
  typename sample_to_stencil_t
  >
struct composeUnsteady<
  ::pressio::rom::lspg::impl::HyperReduced,
  mpl::enable_if_t<
    ::pressio::rom::constraints::continuous_time_system_with_user_provided_apply_jacobian<
      fom_system_type, typename decoder_type::jacobian_type>::value
    >,
  stepper_tag, fom_system_type, decoder_type, lspg_state_type, sample_to_stencil_t>
{
  static_assert(valid_stepper_tag_continuous_time_api<stepper_tag>::value,"");

  // currently only enable when fom types are eigen and pybind11 data structures
  using fom_state_type = typename fom_system_type::state_type;
  static_assert(
#ifdef PRESSIO_ENABLE_TPL_EIGEN
  ::pressio::containers::predicates::is_vector_eigen<fom_state_type>::value
#endif
#if defined(PRESSIO_ENABLE_TPL_EIGEN) and defined(PRESSIO_ENABLE_TPL_PYBIND11)
  or
#endif
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
   ::pressio::containers::predicates::is_array_pybind<fom_state_type>::value
#endif
   , "Unsteady hyper-reduced LSPG with continuous-time API with explicit \
mapping from stencil to sample mesh is currently only enabled for eigen types or pressio4py.");

  using type =
    ::pressio::rom::lspg::impl::unsteady::HyperReducedProblemContinuousTimeApi<
    stepper_tag, fom_system_type, lspg_state_type, decoder_type, sample_to_stencil_t, void>;
};

/***
    preconditioned hyper-reduced lspg with pressio ops and provided mapping from sample to stencil mesh.
    currently only enabled for shared-mem data structures: eigen and pybind11
***/
template<
  typename stepper_tag,
  typename fom_system_type,
  typename decoder_type,
  typename lspg_state_type,
  typename precond_type,
  typename sample_to_stencil_t
  >
struct composeUnsteady<
  ::pressio::rom::lspg::impl::PreconditionedHyperReduced,
  mpl::enable_if_t<
    ::pressio::rom::constraints::continuous_time_system_with_user_provided_apply_jacobian<
      fom_system_type, typename decoder_type::jacobian_type>::value
    >,
  stepper_tag, fom_system_type, decoder_type, lspg_state_type, precond_type, sample_to_stencil_t>
{
  static_assert(valid_stepper_tag_continuous_time_api<stepper_tag>::value,"");

  // currently only enable when fom types are eigen and pybind11 data structures
  using fom_state_type = typename fom_system_type::state_type;
  static_assert(
#ifdef PRESSIO_ENABLE_TPL_EIGEN
  ::pressio::containers::predicates::is_vector_eigen<fom_state_type>::value
#endif
#if defined(PRESSIO_ENABLE_TPL_EIGEN) and defined(PRESSIO_ENABLE_TPL_PYBIND11)
  or
#endif
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
   ::pressio::containers::predicates::is_array_pybind<fom_state_type>::value
#endif
   , "Unsteady preconditioned hyper-reduced LSPG with continuous-time API with explicit \
mapping from stencil to sample mesh is currently only enabled for eigen types or pressio4py.");

  using type =
    ::pressio::rom::lspg::impl::unsteady::PreconditionedHyperReducedProblemContinuousTimeApi<
    stepper_tag, fom_system_type, lspg_state_type, decoder_type, precond_type, sample_to_stencil_t, void>;
};
#endif


#ifdef PRESSIO_ENABLE_TPL_TRILINOS
/***
    hyper-reduced lspg pressio ops with hyp-red enforced in pressio automatically.
    but this is enabled only for Trilinos data types
***/
template<
  typename stepper_tag,
  typename fom_system_type,
  typename decoder_type,
  typename lspg_state_type
  >
struct composeUnsteady<
  ::pressio::rom::lspg::impl::HyperReduced,
  mpl::enable_if_t<
    ::pressio::rom::constraints::continuous_time_system_with_user_provided_apply_jacobian<
      fom_system_type, typename decoder_type::jacobian_type>::value
    >,
  stepper_tag, fom_system_type, decoder_type, lspg_state_type>
{
  static_assert(valid_stepper_tag_continuous_time_api<stepper_tag>::value,"");

  using fom_state_type = typename fom_system_type::state_type;
  static_assert
  (::pressio::containers::predicates::is_vector_tpetra<fom_state_type>::value or
   ::pressio::containers::predicates::is_vector_epetra<fom_state_type>::value or
   ::pressio::containers::predicates::is_vector_tpetra_block<fom_state_type>::value,
   "Unsteady hyper-reduced LSPG with continuous-time API with automatic \
enforcement for hyp-red is currently only supported for Trilinos data types.");

  using type =
    ::pressio::rom::lspg::impl::unsteady::HyperReducedProblemContinuousTimeApi<
    stepper_tag, fom_system_type, lspg_state_type, decoder_type, void, void>;
};

/***
    preconditioned hyper-reduced lspg pressio ops with hyp-red enforced in pressio automatically.
    but this is enabled only for Trilinos data types
***/
template<
  typename stepper_tag,
  typename fom_system_type,
  typename decoder_type,
  typename lspg_state_type,
  typename precond_type
  >
struct composeUnsteady<
  ::pressio::rom::lspg::impl::PreconditionedHyperReduced,
  mpl::enable_if_t<
    ::pressio::rom::constraints::continuous_time_system_with_user_provided_apply_jacobian<
      fom_system_type, typename decoder_type::jacobian_type>::value
    >,
  stepper_tag, fom_system_type, decoder_type, lspg_state_type, precond_type>
{
  static_assert(valid_stepper_tag_continuous_time_api<stepper_tag>::value,"");

  using fom_state_type = typename fom_system_type::state_type;
  static_assert
  (::pressio::containers::predicates::is_vector_tpetra<fom_state_type>::value or
   ::pressio::containers::predicates::is_vector_epetra<fom_state_type>::value or
   ::pressio::containers::predicates::is_vector_tpetra_block<fom_state_type>::value,
   "Unsteady preconditioned hyper-reduced LSPG with continuous-time API with automatic \
enforcement for hyp-red is currently only supported for Trilinos data types.");

  using type =
    ::pressio::rom::lspg::impl::unsteady::PreconditionedHyperReducedProblemContinuousTimeApi<
    stepper_tag, fom_system_type, lspg_state_type, decoder_type, precond_type, void, void>;
};
#endif


///////////////////////////////////////////////////////
//////       DISCRETE TIME API
///////////////////////////////////////////////////////

/***
     default lspg
***/
template<
  typename stepper_tag,
  typename fom_system_type,
  typename decoder_type,
  typename lspg_state_type,
  typename ...Args
  >
struct composeUnsteady<
  ::pressio::rom::lspg::impl::Default,
  mpl::enable_if_t<
    ::pressio::rom::constraints::discrete_time_system_with_user_provided_apply_jacobian<
      fom_system_type, typename decoder_type::jacobian_type>::value
    >,
  stepper_tag, fom_system_type, decoder_type, lspg_state_type, Args...>
{
  static_assert
  (std::is_same< stepper_tag, ::pressio::ode::implicitmethods::Arbitrary>::value,
   "LSPG with discrete-time API currently accepts Arbitrary stepper");

  using type =
    ::pressio::rom::lspg::impl::unsteady::DefaultProblemDiscreteTimeApi<
    stepper_tag, fom_system_type, lspg_state_type, decoder_type, Args...>;
};

/***
     hyper-reduced lspg
     (for now this is just an alias of the default since
     there is no major difference with the default case with discrete-time api
     because the user is supposed to assemble all operators anyway)
***/
template<
  typename stepper_tag,
  typename fom_system_type,
  typename decoder_type,
  typename lspg_state_type,
  typename ...Args
  >
struct composeUnsteady<
  ::pressio::rom::lspg::impl::HyperReduced,
  mpl::enable_if_t<
    ::pressio::rom::constraints::discrete_time_system_with_user_provided_apply_jacobian<
      fom_system_type, typename decoder_type::jacobian_type>::value
    >,
  stepper_tag, fom_system_type, decoder_type, lspg_state_type, Args...
  >
  : composeUnsteady<::pressio::rom::lspg::impl::Default, void,
		     stepper_tag, fom_system_type, decoder_type, lspg_state_type, Args...>
{
  using base_t = composeUnsteady
    <::pressio::rom::lspg::impl::Default, void,
      stepper_tag, fom_system_type, decoder_type, lspg_state_type, Args...>;

  using typename base_t::type;
};

/***
    Preconditioned lspg
***/
template<
  typename stepper_tag,
  typename fom_system_type,
  typename decoder_type,
  typename lspg_state_type,
  typename precond_type,
  typename ...Args
  >
struct composeUnsteady<
  ::pressio::rom::lspg::impl::Preconditioned,
  mpl::enable_if_t<
    ::pressio::rom::constraints::discrete_time_system_with_user_provided_apply_jacobian<
      fom_system_type, typename decoder_type::jacobian_type>::value
    >,
  stepper_tag, fom_system_type, decoder_type, lspg_state_type, precond_type, Args...>
{
  static_assert
  (std::is_same< stepper_tag, ::pressio::ode::implicitmethods::Arbitrary>::value,
   "LSPG with discrete-time API currently accepts Arbitrary stepper");

  using type =
    ::pressio::rom::lspg::impl::unsteady::PreconditionedProblemDiscreteTimeApi<
    stepper_tag, fom_system_type, lspg_state_type, decoder_type, precond_type,Args...>;
};

/***
    masked lspg
***/
template<
  typename stepper_tag,
  typename fom_system_type,
  typename decoder_type,
  typename lspg_state_type,
  typename masker_type,
  typename ...Args
  >
struct composeUnsteady<
  ::pressio::rom::lspg::impl::Masked,
  mpl::enable_if_t<
    ::pressio::rom::constraints::discrete_time_system_with_user_provided_apply_jacobian<
      fom_system_type, typename decoder_type::jacobian_type>::value
    >,
  stepper_tag, fom_system_type, decoder_type, lspg_state_type, masker_type, Args...>
{
  static_assert
  (std::is_same< stepper_tag, ::pressio::ode::implicitmethods::Arbitrary>::value,
   "LSPG with discrete-time API currently accepts Arbitrary stepper");

  using type =
    ::pressio::rom::lspg::impl::unsteady::MaskedProblemDiscreteTimeApi<
    stepper_tag, fom_system_type, lspg_state_type, decoder_type, masker_type, Args...>;
};

}}}}
#endif  // ROM_LSPG_IMPL_UNSTEADY_ROM_COMPOSE_UNSTEADY_LSPG_IMPL_HPP_
