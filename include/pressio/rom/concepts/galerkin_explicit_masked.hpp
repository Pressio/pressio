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

#ifndef ROM_CONSTRAINTS_ROM_UNSTEADY_EXPLICIT_GALERKIN_HYPREDUCEABLE_MASKABLE_WITH_HPP_
#define ROM_CONSTRAINTS_ROM_UNSTEADY_EXPLICIT_GALERKIN_HYPREDUCEABLE_MASKABLE_WITH_HPP_

#include "helpers.hpp"

#ifdef PRESSIO_ENABLE_CXX20

// this is here so that we can clearly show it in the
// doc via rst literal include directive
namespace pressio{ namespace rom{ namespace galerkin{ namespace unsteadyexplicit{

template <
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType,
  class HyperReducerType
  >
concept ComposableIntoHyperReducedMaskedProblem =
     PossiblyAffineTrialColumnSubspace<TrialSubspaceType>
  && SemiDiscreteFom<FomSystemType>
  && std::same_as<
       typename TrialSubspaceType::full_state_type,
       typename FomSystemType::state_type>
  // the masker acts on the fom
  && requires(const MaskerType & masker,
	      const typename FomSystemType::right_hand_side_type & rhsFom)
  {
    /*
      Define:
      - B     : the subspace's basis, i.e. subspace.basisOfTranslatedSpace()
      - rhsFom  : FOM right_hand_side instance
      - rhsFomMasked : the result of the masker acting on rhsFom
    */

     { masker.createApplyMaskResult(rhsFom)  } -> std::copy_constructible;
     { masker(rhsFom,
	      std::declval<decltype(masker.createApplyMaskResult(rhsFom)) &>()
	      ) } -> std::same_as<void>;
  }
  //
  // the hyper-reducer acts on the masked operators
  && requires(
	const HyperReducerType & hyperReducer,
	const concepts::impl::mask_action_t<
	  MaskerType,
	  typename FomSystemType::right_hand_side_type> & rhsFomMasked,
	const typename FomSystemType::time_type & evalTime,
	impl::explicit_galerkin_default_reduced_right_hand_side_t<TrialSubspaceType> & rhsGal)
  {

    /*
      Define:
      - rhsFomMasked : the result of the masker acting on rhsFom
      - rhsGal: the reduced Galerkin right_hand_side
    */

    { hyperReducer(rhsFomMasked, evalTime, rhsGal) } -> std::same_as<void>;
  };

}}}} //end namespace pressio::rom::galerkin::unsteadyexplicit


// #else

// namespace pressio{ namespace rom{ namespace galerkin{ namespace steady{

// template<class T, class SubspaceType, class enable = void>
// struct ProjectableOnPossiblyAffineSubspace : std::false_type{};

// template<class T, class SubspaceType>
// struct ProjectableOnPossiblyAffineSubspace<
//   T, SubspaceType,
//   mpl::enable_if_t<
//        SteadyFomWithJacobianAction<T, SubspaceType>::value
//     && PossiblyAffineTrialColumnSubspace<SubspaceType>::value
//     //
//     && std::is_void<
// 	  decltype
// 	  (
// 	  ::pressio::ops::product
// 	  (::pressio::transpose(),
//            std::declval<typename ::pressio::Traits<typename SubspaceType::basis_matrix_type>::scalar_type const &>(),
// 	   std::declval<typename SubspaceType::basis_matrix_type const &>(),
// 	   std::declval<typename T::residual_type const &>(),
// 	   std::declval<typename Traits<typename SubspaceType::reduced_state_type>::scalar_type const &>(),
// 	   std::declval<typename ::pressio::rom::SteadyGalerkinDefaultOperatorsTraits<
// 	     typename SubspaceType::reduced_state_type>::reduced_residual_type &>()
// 	   )
// 	   )
//     >::value
//     //
//     && std::is_void<
// 	  decltype
// 	  (
// 	  ::pressio::ops::product
// 	  (::pressio::transpose(), ::pressio::nontranspose(),
//            std::declval<typename Traits<typename SubspaceType::basis_matrix_type>::scalar_type const& >(),
// 	   std::declval<typename SubspaceType::basis_matrix_type const &>(),
// 	   std::declval<concepts::impl::fom_jacobian_action_t<T, SubspaceType> const &>(),
// 	   std::declval<typename Traits<typename SubspaceType::reduced_state_type>::scalar_type const&>(),
// 	   std::declval<typename ::pressio::rom::SteadyGalerkinDefaultOperatorsTraits<
//              typename SubspaceType::reduced_state_type>::reduced_jacobian_type &
// 	   >()
// 	   )
// 	   )
//     >::value
//    >
//   > : std::true_type{};

// }}}} //end namespace pressio::rom::galerkin::steady
#endif

#endif  // ROM_CONSTRAINTS_ROM_FOM_SYSTEM_CONTINUOUS_TIME_HPP_
