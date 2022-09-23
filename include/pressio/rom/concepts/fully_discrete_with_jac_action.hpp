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

#ifndef ROM_CONSTRAINTS_ROM_ALL_CONCEPTS_HPP_
#define ROM_CONSTRAINTS_ROM_ALL_CONCEPTS_HPP_

namespace pressio{ namespace rom{

// template<class T, int NumStates, class ManifoldJacType, class enable = void>
// struct FullyDiscreteSystemWithJacobianAction : std::false_type{};

// template<class T, int NumStates, class ManifoldJacType>
// struct FullyDiscreteSystemWithJacobianAction<
//   T, NumStates, ManifoldJacType,
//   mpl::enable_if_t<
//        ::pressio::has_time_typedef<T>::value
//     && ::pressio::has_state_typedef<T>::value
//     && ::pressio::has_discrete_residual_typedef<T>::value
//     //
//     && mpl::is_same<
// 	 typename T::discrete_residual_type,
// 	 decltype(std::declval<T const>().createDiscreteTimeResidual())
// 	 >::value
//     && !std::is_void<
// 	decltype
// 	(
// 	 std::declval<T const>().createResultOfDiscreteTimeJacobianActionOn
// 	 (
// 	  std::declval<ManifoldJacType const &>()
// 	  )
// 	 )
// 	>::value

//     && ::pressio::rom::has_const_discrete_residual_jacobian_action_method<
// 	 T, NumStates,
//          typename ::pressio::ode::StepCount::value_type,
// 	 typename T::time_type,
// 	 typename T::state_type,
// 	 typename T::discrete_residual_type,
// 	 ManifoldJacType,
// 	 decltype
// 	 (
// 	 std::declval<T const>().createResultOfDiscreteTimeJacobianActionOn
// 	 (
// 	 std::declval<ManifoldJacType const &>()
// 	 )
// 	 )
// 	 >::value
//     >
//   > : std::true_type{};

}}
#endif  // ROM_CONSTRAINTS_ROM_FOM_SYSTEM_CONTINUOUS_TIME_HPP_
