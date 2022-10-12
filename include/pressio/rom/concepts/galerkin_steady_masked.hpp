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

#ifndef ROM_CONCEPTS_GALERKIN_STEADY_MASKED_HPP_
#define ROM_CONCEPTS_GALERKIN_STEADY_MASKED_HPP_

#include "helpers.hpp"

namespace pressio{ namespace rom{ namespace galerkin{ namespace steady{

#ifdef PRESSIO_ENABLE_CXX20
template <
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType,
  class HyperReducerType>
concept ComposableIntoHyperReducedMaskedProblem =
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
       MaskerType>
  /*
    requirements on the hyper-reducer:
    must be applicable to the masked operators
  */
  && requires(
	const HyperReducerType & hyperReducer,
	const concepts::impl::mask_action_t<
	  MaskerType,
	  typename FomSystemType::residual_type> & rFomMasked,
	const concepts::impl::mask_action_t<
	  MaskerType,
	  concepts::impl::fom_jacobian_action_on_trial_space_t<
	    FomSystemType, TrialSubspaceType>> & JaFomMasked,
	impl::steady_galerkin_default_reduced_residual_t<TrialSubspaceType> & rGal,
	impl::steady_galerkin_default_reduced_jacobian_t<TrialSubspaceType> & JGal)
  {

    /*
      Define:
      - rFomMasked  : the result of the masker acting on rFom
      - JaFomMasked : the result of the masker acting on JaFom
      - rGal        : the reduced Galerkin residual
      - JGal        : the reduced Galerkin Jacobian
    */

    { hyperReducer(rFomMasked,  rGal) } -> std::same_as<void>;
    { hyperReducer(JaFomMasked, JGal) } -> std::same_as<void>;
  };
#endif // PRESSIO_ENABLE_CXX20

}}}} //end namespace pressio::rom::galerkin::steady


#if not defined PRESSIO_ENABLE_CXX20
namespace pressio{ namespace rom{ namespace galerkin{ namespace steady{

template<
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType,
  class HyperReducerType,
  class enable = void>
struct ComposableIntoHyperReducedMaskedProblem : std::false_type{};

template<
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType,
  class HyperReducerType>
struct ComposableIntoHyperReducedMaskedProblem<
  TrialSubspaceType, FomSystemType, MaskerType, HyperReducerType,
  mpl::enable_if_t<
    PossiblyAffineTrialColumnSubspace<TrialSubspaceType>::value
    && SteadyFomWithJacobianAction<FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value
    //
    && std::is_same<
       typename TrialSubspaceType::full_state_type,
       typename FomSystemType::state_type>::value
    && MaskableWith<
         typename FomSystemType::residual_type,
         MaskerType>::value
    && MaskableWith<
         concepts::impl::fom_jacobian_action_on_trial_space_t<FomSystemType, TrialSubspaceType>,
         MaskerType>::value
    && std::is_void<
       decltype
       (
	std::declval<HyperReducerType const>()
	(
	  std::declval<
	    concepts::impl::mask_action_t<MaskerType, typename FomSystemType::residual_type>
	    const &>(),
	  std::declval<impl::steady_galerkin_default_reduced_residual_t<TrialSubspaceType> & >()
	 )
	)
      >::value
    //
    && std::is_void<
       decltype
       (
	std::declval<HyperReducerType const>()
	(
	  std::declval<
	     concepts::impl::mask_action_t<
	       MaskerType,
	       concepts::impl::fom_jacobian_action_on_trial_space_t<FomSystemType, TrialSubspaceType>>
	     const &>(),
	  std::declval<impl::steady_galerkin_default_reduced_jacobian_t<TrialSubspaceType> &>()
	 )
	)
      >::value
    >
  > : std::true_type{};

}}}} //end namespace pressio::rom::galerkin::steady

#endif

#endif  // ROM_CONCEPTS_GALERKIN_STEADY_MASKED_HPP_
