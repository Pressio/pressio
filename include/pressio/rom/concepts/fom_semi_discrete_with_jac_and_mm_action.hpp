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

#ifndef ROM_CONCEPTS_FOM_SEMI_DISCRETE_WITH_JAC_AND_MM_ACTION_HPP_
#define ROM_CONCEPTS_FOM_SEMI_DISCRETE_WITH_JAC_AND_MM_ACTION_HPP_

#include "helpers.hpp"

namespace pressio{ namespace rom{

#ifdef PRESSIO_ENABLE_CXX20
template<class T, class JacobianActionOperandType, class MassMatrixActionOperandType>
concept SemiDiscreteFomWithJacobianAndMassMatrixAction =
     SemiDiscreteFomWithJacobianAction<T, JacobianActionOperandType>
  && SemiDiscreteFomWithMassMatrixAction<T, MassMatrixActionOperandType>;
#endif // PRESSIO_ENABLE_CXX20

}} // end namespace pressio::rom


#if not defined PRESSIO_ENABLE_CXX20
namespace pressio{ namespace rom{

template<
  class T,
  class JacobianActionOperandType,
  class MassMatrixActionOperandType,
  class enable = void>
struct SemiDiscreteFomWithJacobianAndMassMatrixAction : std::false_type{};

template<class T, class JacobianActionOperandType, class MassMatrixActionOperandType>
struct SemiDiscreteFomWithJacobianAndMassMatrixAction<
  T,  JacobianActionOperandType, MassMatrixActionOperandType,
  mpl::enable_if_t<
     SemiDiscreteFomWithJacobianAction<T, JacobianActionOperandType>::value
  && SemiDiscreteFomWithMassMatrixAction<T, MassMatrixActionOperandType>::value
  >
  > : std::true_type{};

}} // end namespace pressio::rom
#endif

#endif  // ROM_CONCEPTS_FOM_SEMI_DISCRETE_WITH_JAC_ACTION_HPP_