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

#ifndef ROM_CONSTRAINTS_ROM_STEADY_GALERKIN_DEFAULT_PROJECTABLE_WITH_CONCEPTS_HPP_
#define ROM_CONSTRAINTS_ROM_STEADY_GALERKIN_DEFAULT_PROJECTABLE_WITH_CONCEPTS_HPP_

#include "helpers.hpp"

namespace pressio{ namespace rom{
namespace galerkin{ namespace steady{

#ifdef PRESSIO_ENABLE_CXX20

template <class T, class SubspaceType>
concept ProjectableOnPossiblyAffineSubspace =
      SteadyFomWithJacobianAction<T, SubspaceType>
   && PossiblyAffineTrialColumnSubspace<SubspaceType>
   && requires(
	const SubspaceType & subspace,
	typename SteadyGalerkinDefaultOperatorsTraits<
	  typename SubspaceType::reduced_state_type>::reduced_residual_type & rGal,
	typename SteadyGalerkinDefaultOperatorsTraits<
	  typename SubspaceType::reduced_state_type>::reduced_jacobian_type & JGal,
	typename Traits<typename SubspaceType::basis_matrix_type>::scalar_type alpha,
	typename Traits<typename SubspaceType::reduced_state_type>::scalar_type beta)
  {

    /*
      Define:
      - B     : the subspace's basis, i.e. subspace.basisOfTranslatedSpace()
      - rFom  : FOM residual instance
      - JFom  : FOM jacobian instance
      - JaFom : the action of the FOM jacobian on B, i.e. JFom*B
      - alpha : a scalar
      - beta  : a scalar

      - rGal: the reduced Galerkin residual
      - JGal: the reduced Galerkin Jacobian
    */

    /* rGal = beta*rGal + alpha * B^T * rFom */
    pressio::ops::product
       (::pressio::transpose(),
	alpha,
	subspace.basisOfTranslatedSpace(),
	std::declval<const typename T::residual_type &>(),
	beta,
	rGal);

    /* JGal = beta*JGal + alpha * B^T * JaFom * B */
    pressio::ops::product
       (::pressio::transpose(), ::pressio::nontranspose(),
	alpha,
	subspace.basisOfTranslatedSpace(),
	std::declval<
	  const decltype(std::declval<T const &>().createApplyJacobianResult
	      (subspace.basisOfTranslatedSpace())) &
	>(),
	beta,
	JGal);
  };

#else

template<class T, class SubspaceType, class enable = void>
struct ProjectableOnPossiblyAffineSubspace : std::false_type{};

template<class T, class SubspaceType>
struct ProjectableOnPossiblyAffineSubspace<
  T, SubspaceType,
  mpl::enable_if_t<
       SteadyFomWithJacobianAction<T, SubspaceType>::value
    && PossiblyAffineTrialColumnSubspace<SubspaceType>::value
    //
    && std::is_void<
	  decltype
	  (
	  ::pressio::ops::product
	  (::pressio::transpose(),
           std::declval<typename ::pressio::Traits<typename SubspaceType::basis_matrix_type>::scalar_type const &>(),
	   std::declval<typename SubspaceType::basis_matrix_type const &>(),
	   std::declval<typename T::residual_type const &>(),
	   std::declval<typename Traits<typename SubspaceType::reduced_state_type>::scalar_type const &>(),
	   std::declval<typename ::pressio::rom::SteadyGalerkinDefaultOperatorsTraits<
	     typename SubspaceType::reduced_state_type>::reduced_residual_type &>()
	   )
	   )
    >::value
    //
    && std::is_void<
	  decltype
	  (
	  ::pressio::ops::product
	  (::pressio::transpose(), ::pressio::nontranspose(),
           std::declval<typename Traits<typename SubspaceType::basis_matrix_type>::scalar_type const& >(),
	   std::declval<typename SubspaceType::basis_matrix_type const &>(),
	   std::declval<concepts::impl::FomJacActionResult_t<T, SubspaceType> const &>(),
	   std::declval<typename Traits<typename SubspaceType::reduced_state_type>::scalar_type const&>(),
	   std::declval<typename ::pressio::rom::SteadyGalerkinDefaultOperatorsTraits<
             typename SubspaceType::reduced_state_type>::reduced_jacobian_type &
	   >()
	   )
	   )
    >::value
   >
  > : std::true_type{};

#endif

}}}}// end namespace galerkin::steaay
#endif  // ROM_CONSTRAINTS_ROM_FOM_SYSTEM_CONTINUOUS_TIME_HPP_
