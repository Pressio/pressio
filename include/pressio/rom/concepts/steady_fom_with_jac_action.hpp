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

#ifndef ROM_CONSTRAINTS_ROM_STEADY_FOM_WITH_JAC_ACTION_CONCEPT_HPP_
#define ROM_CONSTRAINTS_ROM_STEADY_FOM_WITH_JAC_ACTION_CONCEPT_HPP_

#include "helpers.hpp"

namespace pressio{ namespace rom{

#ifdef PRESSIO_ENABLE_CXX20

template <class T, class TrialSubspaceType>
concept SteadyFomWithJacobianAction =
      PossiblyAffineTrialColumnSubspace<TrialSubspaceType>
   && std::copy_constructible<typename T::state_type>
   && std::copy_constructible<typename T::residual_type>
   && std::same_as<
       typename pressio::Traits<typename T::state_type>::scalar_type,
       typename pressio::Traits<typename T::residual_type>::scalar_type>
   && requires(const T & A,
               const typename T::state_type & state,
               typename T::residual_type & residual,
               const typename TrialSubspaceType::basis_matrix_type & basisMatrix)
   {
     { A.createResidual()          } -> std::same_as<typename T::residual_type>;
     { A.residual(state, residual) } -> std::same_as<void>;

     { A.createApplyJacobianResult(basisMatrix) } -> std::copy_constructible;

     { A.applyJacobian(state, basisMatrix,
                       std::declval<decltype(A.createApplyJacobianResult(basisMatrix)) &>()
                       )} -> std::same_as<void>;
   };


//&& std::same_as<
// typename pressio::Traits<decltype(A.createApplyJacobianResult(basisMatrix))>::scalar_type,
//typename pressio::Traits<typename T::residual_type>::scalar_type>;

#else

template<class T, class TrialSubspaceType, class enable = void>
struct SteadyFomWithJacobianAction : std::false_type{};

template<class T, class TrialSubspaceType>
struct SteadyFomWithJacobianAction<
  T, TrialSubspaceType,
  mpl::enable_if_t<
       PossiblyAffineTrialColumnSubspace<TrialSubspaceType>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_residual_typedef<T>::value
    //
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::residual_type>::value
    && ::pressio::VectorSpaceElementsWithSameField<
      typename T::state_type, typename T::residual_type
      >::value
    //
    && ::pressio::rom::has_const_create_residual_method_return_result<
      T, typename T::residual_type >::value
    && ::pressio::rom::has_const_residual_method_accept_state_result_return_void<
      T, typename T::state_type, typename T::residual_type>::value
    //
    && ::pressio::rom::has_const_create_apply_jacobian_result_method_accept_operand_return_result<
	 T, typename TrialSubspaceType::basis_matrix_type>::value
    && ::pressio::rom::has_const_apply_jacobian_method_accept_state_operand_result_return_void<
	 T, typename T::state_type, typename TrialSubspaceType::basis_matrix_type,
	 concepts::impl::FomJacActionResult_t<T, TrialSubspaceType>
	 >::value
    //
    && std::is_copy_constructible<
	 concepts::impl::FomJacActionResult_t<T, TrialSubspaceType>
	 >::value
    && ::pressio::VectorSpaceElementsWithSameField<
	 typename T::state_type,
  	 concepts::impl::FomJacActionResult_t<T, TrialSubspaceType>
	 >::value
   >
  > : std::true_type{};

#endif

}}
#endif  // ROM_CONSTRAINTS_ROM_FOM_SYSTEM_CONTINUOUS_TIME_HPP_
