/*
//@HEADER
// ************************************************************************
//
// rom_fom_system_continuous_time.hpp
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

#ifndef ROM_CONSTRAINTS_ROM_STEADY_LSPG_MASKABLE_WITH_HPP_
#define ROM_CONSTRAINTS_ROM_STEADY_LSPG_MASKABLE_WITH_HPP_

#include "helpers.hpp"

#ifdef PRESSIO_ENABLE_CXX20

// this is here so that we can clearly show it in the
// doc via rst literal include directive
namespace pressio{ namespace rom{ namespace lspg{ namespace steady{

template <class TrialSubspaceType, class FomSystemType, class MaskerType>
concept ComposableIntoMaskedProblem =
  PossiblyAffineTrialColumnSubspace<TrialSubspaceType>
  && SteadyFomWithJacobianAction<FomSystemType, typename TrialSubspaceType::basis_matrix_type>
  && std::same_as<
       typename TrialSubspaceType::full_state_type,
       typename FomSystemType::state_type>
  /*
    requirements on the masker
    must be applicable to a FOM residual instance as well as
    to the result of the FOM jacobian action on the basis
  */
  && MaskableWith<
       typename FomSystemType::residual_type,
       MaskerType>
  && MaskableWith<
       concepts::impl::fom_jacobian_action_on_trial_space_t<FomSystemType, TrialSubspaceType>,
       MaskerType>;

}}}} //end namespace pressio::rom::lspg::steady




/* leave some white space on purpose so that
   if we make edits above, we don't have to change
   the line numbers included in the rst doc page */

#else

namespace pressio{ namespace rom{ namespace lspg{ namespace steady{

template<
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType,
  class enable = void>
struct ComposableIntoMaskedProblem : std::false_type{};

template<class TrialSubspaceType, class FomSystemType, class MaskerType>
struct ComposableIntoMaskedProblem<
  TrialSubspaceType, FomSystemType, MaskerType,
  mpl::enable_if_t<
    PossiblyAffineTrialColumnSubspace<TrialSubspaceType>::value
    && SteadyFomWithJacobianAction<
        FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value
    && std::is_same<
       typename TrialSubspaceType::full_state_type,
       typename FomSystemType::state_type>::value
    && MaskableWith<
        typename FomSystemType::residual_type, MaskerType>::value
    && MaskableWith<
        concepts::impl::fom_jacobian_action_on_trial_space_t<FomSystemType, TrialSubspaceType>,
        MaskerType>::value
    >
  > : std::true_type{};

}}}} //end namespace pressio::rom::lspg::steady
#endif

#endif  // ROM_CONSTRAINTS_ROM_FOM_SYSTEM_CONTINUOUS_TIME_HPP_
