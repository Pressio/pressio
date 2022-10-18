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

#ifndef ROM_CONCEPTS_LSPG_UNSTEADY_MASKED_HPP_
#define ROM_CONCEPTS_LSPG_UNSTEADY_MASKED_HPP_

#include "helpers.hpp"

namespace pressio{ namespace rom{ namespace lspg{ namespace unsteady{

#ifdef PRESSIO_ENABLE_CXX20
template <class TrialSubspaceType, class FomSystemType, class MaskerType>
concept ComposableIntoMaskedProblem =
     SemiDiscreteFomWithJacobianAction<FomSystemType, typename TrialSubspaceType::basis_matrix_type>
  && ComposableIntoDefaultProblem<TrialSubspaceType, FomSystemType>
  && ::pressio::ops::is_known_data_type<typename FomSystemType::state_type>::value
  && ::pressio::ops::is_known_data_type<typename FomSystemType::right_hand_side_type>::value
  && ::pressio::ops::is_known_data_type<
       concepts::impl::fom_jacobian_action_on_trial_space_t<FomSystemType, TrialSubspaceType>
     >::value
  && MaskableWith<typename FomSystemType::right_hand_side_type, MaskerType>
  && MaskableWith<
       concepts::impl::fom_jacobian_action_on_trial_space_t<FomSystemType, TrialSubspaceType>,
       MaskerType>
  && ::pressio::ops::is_known_data_type<
        concepts::impl::mask_action_t<MaskerType, typename FomSystemType::right_hand_side_type>
     >::value
  && ::pressio::ops::is_known_data_type<
        concepts::impl::mask_action_t<
	  MaskerType,
	  concepts::impl::fom_jacobian_action_on_trial_space_t<FomSystemType, TrialSubspaceType>
	  >
     >::value;
#endif // PRESSIO_ENABLE_CXX20

}}}} //end namespace pressio::rom::lspg::unsteady


#if not defined PRESSIO_ENABLE_CXX20
namespace pressio{ namespace rom{ namespace lspg{ namespace unsteady{

template <class TrialSubspaceType, class FomSystemType, class MaskerType, class = void>
struct ComposableIntoMaskedProblem : std::false_type{};

template <class TrialSubspaceType, class FomSystemType, class MaskerType>
struct ComposableIntoMaskedProblem<
  TrialSubspaceType, FomSystemType, MaskerType,
  mpl::enable_if_t<
    SemiDiscreteFomWithJacobianAction<FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value
    && PossiblyAffineTrialColumnSubspace<TrialSubspaceType>::value
    && std::is_same<
       typename TrialSubspaceType::full_state_type,
      typename FomSystemType::state_type>::value
    && ::pressio::ops::is_known_data_type<typename FomSystemType::state_type>::value
    && ::pressio::ops::is_known_data_type<typename FomSystemType::right_hand_side_type>::value
    && ::pressio::ops::is_known_data_type<
          concepts::impl::fom_jacobian_action_on_trial_space_t<FomSystemType, TrialSubspaceType>
       >::value
    && MaskableWith<typename FomSystemType::right_hand_side_type, MaskerType>::value
    && MaskableWith<
         concepts::impl::fom_jacobian_action_on_trial_space_t<FomSystemType, TrialSubspaceType>,
      MaskerType>::value
    && ::pressio::ops::is_known_data_type<
        concepts::impl::mask_action_t<MaskerType, typename FomSystemType::right_hand_side_type>
       >::value
    && ::pressio::ops::is_known_data_type<
          concepts::impl::mask_action_t<
	    MaskerType,
	    concepts::impl::fom_jacobian_action_on_trial_space_t<FomSystemType, TrialSubspaceType>
	  >
       >::value
    >
  > : std::true_type{};

}}}}
#endif

#endif  // ROM_CONCEPTS_LSPG_UNSTEADY_MASKED_HPP_
