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

#ifndef ROM_CONSTRAINTS_ROM_FULLY_DISCRETE_WITH_JAC_ACTION_CONCEPT_HPP_
#define ROM_CONSTRAINTS_ROM_FULLY_DISCRETE_WITH_JAC_ACTION_CONCEPT_HPP_

namespace pressio{ namespace rom{

#ifdef PRESSIO_ENABLE_CXX20

template<class T, int NumStates, class TrialSubspaceType>
concept FullyDiscreteSystemWithJacobianAction =
      PossiblyAffineTrialColumnSubspace<TrialSubspaceType>
   && (NumStates == 2 || NumStates == 3)
   && std::regular<typename T::time_type>
   && std::totally_ordered<typename T::time_type>
   && std::copy_constructible<typename T::state_type>
   && std::copy_constructible<typename T::discrete_residual_type>
   && std::same_as<
       typename pressio::Traits<typename T::state_type>::scalar_type,
       typename pressio::Traits<typename T::discrete_residual_type>::scalar_type>
   && std::convertible_to<
       typename T::time_type,
       typename pressio::Traits<typename T::state_type>::scalar_type>
  && requires(const T & A,
	      const typename T::state_type & state,
	      const typename TrialSubspaceType::basis_matrix_type & basisMatrix)
   {
     { A.createDiscreteTimeResidual() } -> std::same_as<typename T::discrete_residual_type>;

     { A.createResultOfDiscreteTimeJacobianActionOn(basisMatrix) } -> std::copy_constructible;
   }
   && ::pressio::rom::has_const_discrete_residual_jacobian_action_method<
	T, NumStates,
	typename ::pressio::ode::StepCount::value_type,
	typename T::time_type,
	typename T::state_type,
	typename T::discrete_residual_type,
	typename TrialSubspaceType::basis_matrix_type,
	decltype
	(
	 std::declval<T const>().createResultOfDiscreteTimeJacobianActionOn
	  (
	   std::declval<typename TrialSubspaceType::basis_matrix_type const &>()
	  )
	 )
     >::value;

#else

template<class T, int NumStates, class TrialSubspaceType, class enable = void>
struct FullyDiscreteSystemWithJacobianAction : std::false_type{};

template<class T, int NumStates, class TrialSubspaceType>
struct FullyDiscreteSystemWithJacobianAction<
  T, NumStates, TrialSubspaceType,
  mpl::enable_if_t<
       ::pressio::has_time_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_discrete_residual_typedef<T>::value
    //
    && mpl::is_same<
	 typename T::discrete_residual_type,
	 decltype(std::declval<T const>().createDiscreteTimeResidual())
	 >::value
    //
    && std::is_copy_constructible<
	decltype
	(
	 std::declval<T const>().createResultOfDiscreteTimeJacobianActionOn
	  (
	   std::declval<typename TrialSubspaceType::basis_matrix_type const &>()
	  )
	 )
	>::value

    && ::pressio::rom::has_const_discrete_residual_jacobian_action_method<
	 T, NumStates,
         typename ::pressio::ode::StepCount::value_type,
	 typename T::time_type,
	 typename T::state_type,
	 typename T::discrete_residual_type,
	 typename TrialSubspaceType::basis_matrix_type,
	 decltype
	 (
	  std::declval<T const>().createResultOfDiscreteTimeJacobianActionOn
	   (
	    std::declval<typename TrialSubspaceType::basis_matrix_type const &>()
	   )
	  )
	 >::value
    >
  > : std::true_type{};

#endif

}}
#endif  // ROM_CONSTRAINTS_ROM_FOM_SYSTEM_CONTINUOUS_TIME_HPP_
