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

#ifndef ROM_CONSTRAINTS_ROM_MASKER_HPP_
#define ROM_CONSTRAINTS_ROM_MASKER_HPP_

namespace pressio{ namespace rom{

template<class T, class enable = void>
struct TimeInvariantMasker : std::false_type{};

template<class T>
struct TimeInvariantMasker<
  T,
  mpl::enable_if_t<
       ::pressio::has_operand_typedef<T>::value
    && ::pressio::has_result_typedef<T>::value
    && std::is_same<
	decltype
	(
	 std::declval<T const>().createApplyMaskResult
	 (
	  std::declval<typename T::operand_type const &>()
	  )
	 ),
	typename T::result_type
       >::value
    //
    && std::is_void<
	decltype
	(
	 std::declval<T const>()
	 (
	  std::declval<typename T::operand_type const &>(),
	  std::declval<typename T::result_type &>()
	  )
	 )
	>::value
   >
  > : std::true_type{};

template<
  class T,
  class ReducedResidualType,
  class ReducedJacActionType,
  class enable = void>
struct SteadyGalerkinHyperReductionOperator : std::false_type{};

template<
  class T,
  class ReducedResidualType,
  class ReducedJacActionType>
struct SteadyGalerkinHyperReductionOperator<
  T, ReducedResidualType, ReducedJacActionType,
  mpl::enable_if_t<
       ::pressio::has_residual_operand_typedef<T>::value
    && ::pressio::has_jacobian_action_operand_typedef<T>::value
    && std::is_void<
	decltype
	(
	 std::declval<T const>()
	 (
	  std::declval<typename T::residual_operand_type const &>(),
	  std::declval<ReducedResidualType &>()
	  )
	 )
	>::value
    && std::is_void<
	decltype
	(
	 std::declval<T const>()
	 (
	  std::declval<typename T::jacobian_action_operand_type const &>(),
	  std::declval<ReducedJacActionType &>()
	  )
	 )
	>::value
   >
  > : std::true_type{};

template<
  class T,
  class ReducedResidualType,
  class enable = void>
struct UnsteadyExplicitGalerkinHyperReductionOperator : std::false_type{};

template<
  class T,
  class ReducedResidualType>
struct UnsteadyExplicitGalerkinHyperReductionOperator<
  T, ReducedResidualType,
  mpl::enable_if_t<
       ::pressio::has_time_typedef<T>::value
    && ::pressio::has_operand_typedef<T>::value
    && std::is_void<
	decltype
	(
	 std::declval<T const>()
	 (
	  std::declval<typename T::operand_type const &>(),
	  std::declval<typename T::time_type>(),
	  std::declval<ReducedResidualType &>()
	  )
	 )
	>::value
   >
  > : std::true_type{};

}}
#endif  // ROM_CONSTRAINTS_ROM_FOM_SYSTEM_CONTINUOUS_TIME_HPP_
