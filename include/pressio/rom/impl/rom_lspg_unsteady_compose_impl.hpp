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

// #include "./rom_lspg_decorator_preconditioner.hpp"
// #include "./rom_lspg_decorator_masked.hpp"
#include "./rom_problem_members_common_mixins.hpp"
#include "./rom_lspg_problem_members.hpp"
#include "./rom_lspg_unsteady_policy_residual.hpp"
#include "./rom_lspg_unsteady_policy_jacobian.hpp"
#include "./rom_lspg_unsteady_traits.hpp"
#include "./rom_lspg_unsteady_problem.hpp"

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{

template <typename tag>
struct valid_stepper_tag_continuous_time_api
{
  static_assert
  (std::is_same<tag, ::pressio::ode::BDF1>::value or
   std::is_same<tag, ::pressio::ode::BDF2>::value or
   std::is_same<tag, ::pressio::ode::CrankNicolson>::value,
   "Invalid stepper tag for LSPG unsteady problem with continuous-time API: \
this can be because you used a wrong one, or the current LSPG implementation does \
not support it, or because you added a new ode scheme in the ode package \
but forgot to update the list of implicit tags supported by LSPG which \
currently contains: BDF1, BDF2 or CrankNicolson");

  static constexpr auto value = true;
};

template<
  int id,
  class StepperTag,
  class FomSystemType,
  class DecoderType,
  class LspgStateType,
  class FomReferenceState,
  class ...Args
  >
struct ComposerContTime
{

  static_assert(valid_stepper_tag_continuous_time_api<StepperTag>::value,"");

  static_assert
  (::pressio::rom::decoder<mpl::remove_cvref_t<DecoderType>, LspgStateType>::value,
   "Invalid decoder detected");
  using decoder_jacobian_type = typename mpl::remove_cvref_t<DecoderType>::jacobian_type;

  static_assert
  (::pressio::rom::continuous_time_fom_system_with_user_provided_apply_jacobian<
   mpl::remove_cvref_t<FomSystemType>, decoder_jacobian_type>::value,
   "The FOM system does not satisfy the steady concept");

  static_assert
  (std::is_same<mpl::remove_cvref_t<FomReferenceState>,
   typename mpl::remove_cvref_t<DecoderType>::fom_state_type>::value,
   "The type deduced for the FOM nominal state passed to the create function is not \
  compatible with the FOM state type detected from adapter class");

  using type = ::pressio::rom::lspg::impl::UnsteadyProblem<
    id, StepperTag, FomSystemType, LspgStateType, DecoderType, Args...>;
};

// default
template<class ...Args>
using ComposeDefaultProblemUnsteadyContTime = ComposerContTime<0, Args...>;

// // preconditioned default
// template<class ...Args>
// using ComposePrecDefaultProblemSteady = Composer<2, Args...>;

// // hyperreduced (note that impl-wise this is same as default)
// template<class ...Args>
// using ComposeHyperreducedProblemSteady = ComposeDefaultProblemSteady<Args...>;

// // preconditioned hypred
// template<class ...Args>
// using ComposePrecHypredProblemSteady = ComposePrecDefaultProblemSteady<Args...>;

// // masked
// template<class ...Args>
// using ComposeMaskedProblemSteady = Composer<1, Args...>;

// // preconditioned masked
// template<class ...Args>
// using ComposePrecMaskedProblemSteady = Composer<3, Args...>;

}}}}
#endif  // ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_ROM_COMPOSE_IMPL_HPP_
