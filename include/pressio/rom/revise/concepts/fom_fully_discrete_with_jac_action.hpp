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

#ifndef ROM_CONCEPTS_FOM_FULLY_DISCRETE_WITH_JAC_ACTION_HPP_
#define ROM_CONCEPTS_FOM_FULLY_DISCRETE_WITH_JAC_ACTION_HPP_

#include "helpers.hpp"

namespace pressio{ namespace rom{

#ifdef PRESSIO_ENABLE_CXX20
template<class T, int TotalNumStates, class JacobianActionOperandType>
concept FullyDiscreteSystemWithJacobianAction =
 /*
   allowed number of state
 */
  (TotalNumStates == 2 || TotalNumStates == 3)
  /*
    required nested aliases
  */
  && requires(){
    typename T::time_type;
    typename T::state_type;
    typename T::discrete_residual_type;
  }
  /*
    requirements on the nested aliases
  */
  && std::regular<typename T::time_type>
  && std::totally_ordered<typename T::time_type>
  && std::copy_constructible<typename T::state_type>
  && std::copy_constructible<typename T::discrete_residual_type>
  && all_have_traits_and_same_scalar<
    typename T::state_type, typename T::discrete_residual_type>::value
  && Traits<typename T::state_type>::rank == 1
  && Traits<typename T::discrete_residual_type>::rank == 1
  && std::convertible_to<
    typename T::time_type, scalar_trait_t<typename T::state_type>>
  /*
   requirements on the "create" methods
  */
  && requires(const T & A,
	      const typename T::state_type & state,
	      const JacobianActionOperandType & operand)
  {
    { A.createDiscreteTimeResidual() } -> std::same_as<typename T::discrete_residual_type>;
    { A.createResultOfDiscreteTimeJacobianActionOn(operand) } -> std::copy_constructible;
  }
  && all_have_traits_and_same_scalar<
    typename T::state_type, JacobianActionOperandType,
    concepts::impl::fully_discrete_fom_jacobian_action_t<T, JacobianActionOperandType>
    >::value
  && ::pressio::SameRankAs<
       concepts::impl::fully_discrete_fom_jacobian_action_t<T, JacobianActionOperandType>,
       JacobianActionOperandType>
  /*
    requirements on "evaluation" methods (fix syntax)
    intentionally not lumped with the above for these reasons:
    1. makes sense to split them, since the following depends on the above
    2. helps the compiler with early failure detection
  */
  && ::pressio::rom::has_const_discrete_residual_jacobian_action_method<
      T, TotalNumStates,
      typename ::pressio::ode::StepCount::value_type,
      typename T::time_type,
      typename T::state_type,
      typename T::discrete_residual_type,
      JacobianActionOperandType,
      concepts::impl::fully_discrete_fom_jacobian_action_t<T, JacobianActionOperandType>
    >::value;
#endif // PRESSIO_ENABLE_CXX20

}} // end namespace pressio::rom


#if not defined PRESSIO_ENABLE_CXX20

namespace pressio{ namespace rom{

template<class T, int TotalNumStates, class JacobianActionOperandType, class enable = void>
struct FullyDiscreteSystemWithJacobianAction : std::false_type{};

template<class T, int TotalNumStates, class JacobianActionOperandType>
struct FullyDiscreteSystemWithJacobianAction<
  T, TotalNumStates, JacobianActionOperandType,
  mpl::enable_if_t<
    (TotalNumStates == 2 || TotalNumStates == 3)
    //
    && ::pressio::has_time_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_discrete_residual_typedef<T>::value
    /*
      requirements on the nested aliases
    */
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::discrete_residual_type>::value
    && all_have_traits_and_same_scalar<
        typename T::state_type, typename T::discrete_residual_type>::value
    && Traits<typename T::state_type>::rank == 1
    && Traits<typename T::discrete_residual_type>::rank == 1
    && std::is_convertible<
      typename T::time_type, scalar_trait_t<typename T::state_type>>::value
    /*
      requirements on the "create" methods
    */
    && mpl::is_same<
	 typename T::discrete_residual_type,
	 decltype(std::declval<T const>().createDiscreteTimeResidual())
	 >::value
    && std::is_copy_constructible<
	decltype
	(
	 std::declval<T const>().createResultOfDiscreteTimeJacobianActionOn
	  (
	   std::declval<JacobianActionOperandType const &>()
	  )
	 )
	>::value
    && all_have_traits_and_same_scalar<
         typename T::state_type, JacobianActionOperandType,
         concepts::impl::fully_discrete_fom_jacobian_action_t<T, JacobianActionOperandType>
      >::value
    && all_have_same_rank<
	 concepts::impl::fully_discrete_fom_jacobian_action_t<T, JacobianActionOperandType>,
         JacobianActionOperandType>::value
    /*
      requirements on "evaluation" methods
    */
    && ::pressio::rom::has_const_discrete_residual_jacobian_action_method<
      T, TotalNumStates,
      typename ::pressio::ode::StepCount::value_type,
      typename T::time_type,
      typename T::state_type,
      typename T::discrete_residual_type,
      JacobianActionOperandType,
      concepts::impl::fully_discrete_fom_jacobian_action_t<T, JacobianActionOperandType>
      >::value
    >
  > : std::true_type{};

}}

#endif

#endif  // ROM_CONCEPTS_FOM_FULLY_DISCRETE_WITH_JAC_ACTION_HPP_
