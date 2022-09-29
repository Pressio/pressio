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

#ifdef PRESSIO_ENABLE_CXX20

// this is here so that we can clearly show it in the
// doc via rst literal include directive
namespace pressio{ namespace rom{

template <class T, class JacobianActionOperandType>
concept SteadyFomWithJacobianAction =
  /* must have nested aliases */
  requires(){
    typename T::state_type;
    typename T::residual_type;
  }
  /*
    requirements on the nested aliases
  */
  && std::copy_constructible<typename T::state_type>
  && std::copy_constructible<typename T::residual_type>
  && all_have_traits_and_same_scalar<
       typename T::state_type, typename T::residual_type>::value
  && Traits<typename T::state_type>::rank == 1
  && Traits<typename T::residual_type>::rank == 1
  /*
    requirements on the "create" methods
  */
  && requires(const T & A,
	      const typename T::state_type   & state,
	      const JacobianActionOperandType & operand)
  {
    { A.createResidual() } -> std::same_as<typename T::residual_type>;
    { A.createApplyJacobianResult(operand) } -> std::copy_constructible;
  }
  && all_have_traits_and_same_scalar<
    typename T::state_type, JacobianActionOperandType,
    concepts::impl::fom_jacobian_action_t<T, JacobianActionOperandType>
    >::value
  && ::pressio::SameRankAs<
       concepts::impl::fom_jacobian_action_t<T, JacobianActionOperandType>,
       JacobianActionOperandType>
  /*
    compound requirements on the "evaluation" methods
    we intentionally do not lump this together with the one
    above for these reasons:
    1. it makes sense logically to split them, since the
       requirements of "A.applyJacobian" and "A.residual"
       depend on successfully checking the corresponding create methods
    2. helps the compiler with early failure detection
  */
  && requires(const T & A,
	      const typename T::state_type & state,
	      const JacobianActionOperandType & operand,
	      typename T::residual_type & residual,
	      concepts::impl::fom_jacobian_action_t<T, JacobianActionOperandType> & jAction)
  {
    { A.residual(state, residual) } -> std::same_as<void>;
    { A.applyJacobian(state, operand, jAction) } -> std::same_as<void>;
  };

}} // end namespace pressio::rom









/* leave some white space on purpose so that
   if we make edits above, we don't have to change
   the line numbers included in the rst doc page */

#else

namespace pressio{ namespace rom{

template<class T, class JacobianActionOperandType, class enable = void>
struct SteadyFomWithJacobianAction : std::false_type{};

template<class T, class JacobianActionOperandType>
struct SteadyFomWithJacobianAction<
  T, JacobianActionOperandType,
  mpl::enable_if_t<
       ::pressio::has_state_typedef<T>::value
    && ::pressio::has_residual_typedef<T>::value
    //
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::residual_type>::value
    && all_have_traits_and_same_scalar<
	 typename T::state_type, typename T::residual_type>::value
    && Traits<typename T::state_type>::rank == 1
    && Traits<typename T::residual_type>::rank == 1
    /*
      requirements on the "create" methods
    */
    && ::pressio::rom::has_const_create_residual_method_return_result<
	 T, typename T::residual_type >::value
    && ::pressio::rom::has_const_create_apply_jacobian_result_method_accept_operand_return_result<
	 T, JacobianActionOperandType>::value
    && std::is_copy_constructible<
	 concepts::impl::fom_jacobian_action_t<T, JacobianActionOperandType>
	 >::value
    && all_have_traits_and_same_scalar<
	 typename T::state_type, JacobianActionOperandType,
	 concepts::impl::fom_jacobian_action_t<T, JacobianActionOperandType>
	 >::value
    && all_have_same_rank<
	 concepts::impl::fom_jacobian_action_t<T, JacobianActionOperandType>,
	 JacobianActionOperandType>::value
    /*
      requirements on "evaluation" methods
    */
    && ::pressio::rom::has_const_residual_method_accept_state_result_return_void<
      T, typename T::state_type, typename T::residual_type>::value
    && ::pressio::rom::has_const_apply_jacobian_method_accept_state_operand_result_return_void<
         T, typename T::state_type, JacobianActionOperandType,
	 concepts::impl::fom_jacobian_action_t<T, JacobianActionOperandType>
	 >::value
   >
  > : std::true_type{};

}}

#endif

#endif  // ROM_CONSTRAINTS_ROM_FOM_SYSTEM_CONTINUOUS_TIME_HPP_