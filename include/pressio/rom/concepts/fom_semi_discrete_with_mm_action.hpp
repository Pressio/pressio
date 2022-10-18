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

#ifndef ROM_CONCEPTS_FOM_SEMI_DISCRETE_WITH_MM_ACTION_HPP_
#define ROM_CONCEPTS_FOM_SEMI_DISCRETE_WITH_MM_ACTION_HPP_

#include "helpers.hpp"

namespace pressio{ namespace rom{

#ifdef PRESSIO_ENABLE_CXX20
template <class T, class MassMatrixActionOperandType>
concept SemiDiscreteFomWithMassMatrixAction =
  /* subsumed concept */
  SemiDiscreteFom<T>
  /*
    compound requirements on the "create" method
  */
  && requires(const T & A,
	      const typename T::state_type    & state,
	      const MassMatrixActionOperandType & operand)
  {
    { A.createResultOfMassMatrixActionOn(operand) } -> std::copy_constructible;
  }
  && all_have_traits_and_same_scalar<
      typename T::state_type,
      concepts::impl::fom_mass_matrix_action_t<T, MassMatrixActionOperandType>,
      MassMatrixActionOperandType
     >::value
  && ::pressio::SameRankAs<
       concepts::impl::fom_mass_matrix_action_t<T, MassMatrixActionOperandType>,
       MassMatrixActionOperandType>
  /*
    compound requirements on the "evaluation" method:
    intentionally not lumped with the above for these reasons:
    1. makes sense to split them, since the following depends on the above
    2. helps the compiler with early failure detection
  */
  && requires(const T & A,
	      const typename T::state_type    & state,
	      const typename T::time_type     & evalTime,
	      const MassMatrixActionOperandType & operand,
	      concepts::impl::fom_mass_matrix_action_t<T, MassMatrixActionOperandType> & result)
  {
    { A.applyMassMatrix(state, operand, evalTime, result) } -> std::same_as<void>;
  };
#endif // PRESSIO_ENABLE_CXX20

}} // end namespace pressio::rom


#if not defined PRESSIO_ENABLE_CXX20

namespace pressio{ namespace rom{

template<class T, class MassMatrixActionOperandType, class enable = void>
struct SemiDiscreteFomWithMassMatrixAction : std::false_type{};

template<class T, class MassMatrixActionOperandType>
struct SemiDiscreteFomWithMassMatrixAction<
  T, MassMatrixActionOperandType,
  mpl::enable_if_t<
    SemiDiscreteFom<T>::value
    //
    && std::is_copy_constructible<
      decltype
      (
       std::declval<T const>().createResultOfMassMatrixActionOn
       (
	std::declval<MassMatrixActionOperandType const &>()
	)
       )
      >::value
    //
    && all_have_traits_and_same_scalar<
	 typename T::state_type,
	 concepts::impl::fom_mass_matrix_action_t<T, MassMatrixActionOperandType>,
	 MassMatrixActionOperandType
       >::value
    && all_have_same_rank<
       concepts::impl::fom_mass_matrix_action_t<T, MassMatrixActionOperandType>,
       MassMatrixActionOperandType>::value
    //
    && std::is_void<
       decltype
       (
	std::declval<T const>().applyMassMatrix
	(
	 std::declval<typename T::state_type const&>(),
	 std::declval<MassMatrixActionOperandType const&>(),
	 std::declval<typename T::time_type const &>(),
	 std::declval<concepts::impl::fom_mass_matrix_action_t<T,  MassMatrixActionOperandType> &>()
	 )
	)
       >::value
   >
  > : std::true_type{};

}} // end namespace pressio::rom

#endif

#endif  // ROM_CONCEPTS_FOM_SEMI_DISCRETE_WITH_MM_ACTION_HPP_
