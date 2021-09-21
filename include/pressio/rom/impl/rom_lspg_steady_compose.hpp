/*
//@HEADER
// ************************************************************************
//
// rom_lspg_steady_compose.hpp
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

#ifndef ROM_IMPL_ROM_LSPG_STEADY_COMPOSE_HPP_
#define ROM_IMPL_ROM_LSPG_STEADY_COMPOSE_HPP_

#include "./rom_lspg_steady_decorators.hpp"
#include "./rom_problem_members_common_mixins.hpp"
#include "./rom_lspg_problem_members.hpp"
#include "./rom_lspg_steady_policies.hpp"
#include "./rom_lspg_steady_system.hpp"
#include "./rom_lspg_steady_traits.hpp"
#include "./rom_lspg_steady_problem.hpp"

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{

template<
  int id,
  class FomSystemType,
  class DecoderType,
  class LspgStateType,
  class FomReferenceState,
  class ...Args
  >
struct Composer
{

  static_assert
  (::pressio::rom::decoder<mpl::remove_cvref_t<DecoderType>, LspgStateType>::value,
   "Invalid decoder detected");
  using decoder_jacobian_type = typename mpl::remove_cvref_t<DecoderType>::jacobian_type;

  static_assert
  (::pressio::rom::steady_fom_system_with_user_provided_apply_jacobian<
   mpl::remove_cvref_t<FomSystemType>, decoder_jacobian_type>::value,
   "The FOM system does not satisfy the steady concept");

  static_assert
    (std::is_same<mpl::remove_cvref_t<FomReferenceState>,
     typename mpl::remove_cvref_t<DecoderType>::fom_state_type>::value,
   "The type deduced for the FOM nominal state passed to the create function is not \
  compatible with the FOM state type detected from adapter class");

  using type = ::pressio::rom::lspg::impl::SteadyProblem<
    id, FomSystemType, LspgStateType, DecoderType, Args...>;
};

// default
template<class ...Args>
using ComposeDefaultProblemSteady = Composer<0, Args...>;

// preconditioned default
template<class ...Args>
using ComposePrecDefaultProblemSteady = Composer<2, Args...>;

// hyperreduced (note that impl-wise this is same as default)
template<class ...Args>
using ComposeHyperreducedProblemSteady = ComposeDefaultProblemSteady<Args...>;

// preconditioned hypred
template<class ...Args>
using ComposePrecHypredProblemSteady = ComposePrecDefaultProblemSteady<Args...>;

// masked
template<class ...Args>
using ComposeMaskedProblemSteady = Composer<1, Args...>;

// preconditioned masked
template<class ...Args>
using ComposePrecMaskedProblemSteady = Composer<3, Args...>;

}}}}
#endif  // ROM_IMPL_ROM_LSPG_STEADY_COMPOSE_HPP_
