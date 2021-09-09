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

#ifndef ROM_LSPG_IMPL_ROM_COMPOSE_UNSTEADY_IMPL_HPP_
#define ROM_LSPG_IMPL_ROM_COMPOSE_UNSTEADY_IMPL_HPP_

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "./rom_lspg_unsteady_hypred_updater_trilinos.hpp"
#endif
#include "./rom_lspg_unsteady_discrete_time_decorators.hpp"
#include "./rom_lspg_unsteady_cont_time_decorators.hpp"
#include "./rom_lspg_unsteady_discrete_time_default_system.hpp"
#include "./rom_problem_members_common_mixins.hpp"
#include "./rom_lspg_problem_members.hpp"
#include "./rom_lspg_unsteady_hypred_policy_residual.hpp"
#include "./rom_lspg_unsteady_hypred_policy_jacobian.hpp"
#include "./rom_lspg_unsteady_policy_residual.hpp"
#include "./rom_lspg_unsteady_policy_jacobian.hpp"
#include "./rom_lspg_unsteady_traits.hpp"
#include "./rom_lspg_unsteady_problem.hpp"

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{

template<
  int id,
  class FomSystemType,
  class DecoderType,
  class LspgStateType,
  class FomReferenceState,
  class ...Args
  >
struct ComposerContTime
{

  static_assert
  (::pressio::rom::decoder<mpl::remove_cvref_t<DecoderType>, LspgStateType>::value,
   "Invalid decoder detected");
  using decoder_jacobian_type = typename mpl::remove_cvref_t<DecoderType>::jacobian_type;

  static_assert
  (::pressio::rom::continuous_time_fom_system_with_user_provided_apply_jacobian<
   mpl::remove_cvref_t<FomSystemType>, decoder_jacobian_type>::value,
   "The FOM system does not satisfy the unsteady concept");

  static_assert
  (std::is_same<mpl::remove_cvref_t<FomReferenceState>,
   typename mpl::remove_cvref_t<DecoderType>::fom_state_type>::value,
   "The type deduced for the FOM nominal state passed to the create function is not \
  compatible with the FOM state type detected from adapter class");

  using type = ::pressio::rom::lspg::impl::UnsteadyProblem<
    id, FomSystemType, LspgStateType, DecoderType, Args...>;
};

template<
  int id,
  std::size_t num_states,
  class FomSystemType,
  class DecoderType,
  class LspgStateType,
  class FomReferenceState,
  class ...Args
  >
struct ComposerDiscTime
{

  static_assert
  (::pressio::rom::decoder<mpl::remove_cvref_t<DecoderType>, LspgStateType>::value,
   "Invalid decoder detected");
  using decoder_jacobian_type = typename mpl::remove_cvref_t<DecoderType>::jacobian_type;

  static_assert
  (::pressio::rom::discrete_time_fom_system_with_user_provided_apply_jacobian<
   mpl::remove_cvref_t<FomSystemType>, num_states, decoder_jacobian_type>::value,
   "The FOM system does not satisfy the discrete-time concept");

  static_assert
  (std::is_same<mpl::remove_cvref_t<FomReferenceState>,
   typename mpl::remove_cvref_t<DecoderType>::fom_state_type>::value,
   "The type deduced for the FOM nominal state passed to the create function is not \
  compatible with the FOM state type detected from adapter class");

  using type = ::pressio::rom::lspg::impl::UnsteadyProblem<
    id, FomSystemType, LspgStateType, DecoderType,
    ::pressio::ode::StepperTotalNumberOfStates<num_states>,  Args...>;
};

// default
template<class ...Args>
using ComposeDefaultProblemContTime = ComposerContTime<0, Args...>;

template<std::size_t n, class ...Args>
using ComposeDefaultProblemDiscTime = ComposerDiscTime<1, n, Args...>;

// preconditioned default
template<class ...Args>
using ComposePrecDefaultProblemContTime = ComposerContTime<2, Args...>;

template<std::size_t n, class ...Args>
using ComposePrecDefaultProblemDiscTime = ComposerDiscTime<3, n, Args...>;

// masked
template<class ...Args>
using ComposeMaskedProblemContTime = ComposerContTime<4, Args...>;

template<std::size_t n, class ...Args>
using ComposeMaskedProblemDiscTime = ComposerDiscTime<5, n, Args...>;

// preconditioned masked
template<class ...Args>
using ComposePrecMaskedProblemContTime = ComposerContTime<6, Args...>;

template<std::size_t n, class ...Args>
using ComposePrecMaskedProblemDiscTime = ComposerDiscTime<7, n, Args...>;

// hyper-reduced
template<class ...Args>
using ComposeHypRedProblemContTime = ComposerContTime<8, Args...>;

// for discrete-time hyperreduced, impl is same as default
template<std::size_t n, class ...Args>
using ComposeHypRedProblemDiscTime = ComposerDiscTime<1, n, Args...>;


// #ifdef PRESSIO_ENABLE_TPL_TRILINOS
// template<>
// struct SelectHypRedUpdaterTrilinos
// {
// }
// #endif


}}}}
#endif  // ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_ROM_COMPOSE_IMPL_HPP_
