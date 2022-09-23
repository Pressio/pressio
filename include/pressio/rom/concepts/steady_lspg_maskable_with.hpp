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

namespace pressio{ namespace rom{
namespace lspg{ namespace steady{

#ifdef PRESSIO_ENABLE_CXX20

template <class T, class MaskerType, class TrialSubspaceType>
concept MaskableWith =
      PossiblyAffineTrialColumnSubspace<TrialSubspaceType>
   && SteadyFomWithJacobianAction<T, TrialSubspaceType>
   && requires(const typename TrialSubspaceType::basis_matrix_type & basisMatrix,
               const typename T::residual_type & fomRe,
               const decltype(std::declval<T const &>().createApplyJacobianResult(basisMatrix)) & fomJa,
               const MaskerType & masker)
  {
    // a masker must be able to operate on the fom
    { masker.createApplyMaskResult(fomRe) } -> std::copy_constructible;
    { masker.createApplyMaskResult(fomJa) } -> std::copy_constructible;
    { masker(fomRe,
	     std::declval<decltype(masker.createApplyMaskResult(fomRe)) &>()
	     ) } -> std::same_as<void>;
    { masker(fomJa,
	     std::declval<decltype(masker.createApplyMaskResult(fomJa)) &>()
	     ) } -> std::same_as<void>;
  };

#else

template<class T, class MaskerType, class TrialSubspaceType, class enable = void>
struct MaskableWith : std::false_type{};

template<class T, class MaskerType, class TrialSubspaceType>
struct MaskableWith<
  T, MaskerType, TrialSubspaceType,
  mpl::enable_if_t<
       PossiblyAffineTrialColumnSubspace<TrialSubspaceType>::value
    && SteadyFomWithJacobianAction<T, TrialSubspaceType>::value
    //
    && !std::is_void<
	 concepts::impl::CreateApplyMaskResult_t<MaskerType, typename T::residual_type>
	 >::value
    && !std::is_void<
	 concepts::impl::CreateApplyMaskResult_t<
	   MaskerType, concepts::impl::FomJacActionResult_t<T, TrialSubspaceType>>
	 >::value
    //
    &&  std::is_void<
	decltype
	(
	 std::declval<MaskerType const>()
	 (
	   std::declval<typename T::residual_type const &>(),
	   std::declval<concepts::impl::CreateApplyMaskResult_t<MaskerType, typename T::residual_type> &>()
	  )
	 )
	 >::value
    &&  std::is_void<
	decltype
	(
	 std::declval<MaskerType const>()
	 (
	   std::declval<concepts::impl::FomJacActionResult_t<T, TrialSubspaceType> const &>(),
	   std::declval<concepts::impl::CreateApplyMaskResult_t<
	     MaskerType, concepts::impl::FomJacActionResult_t<T, TrialSubspaceType>> &>()
	  )
	 )
        >::value
    >
  > : std::true_type{};

#endif

}}}}// end namespace
#endif  // ROM_CONSTRAINTS_ROM_FOM_SYSTEM_CONTINUOUS_TIME_HPP_
