/*
//@HEADER
// ************************************************************************
//
// rom_decoder.hpp
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

#ifndef ROM_CONSTRAINTS_ROM_DECODER_HPP_
#define ROM_CONSTRAINTS_ROM_DECODER_HPP_

namespace pressio{ namespace rom{

/*
 * T is a legitimate decoder if:
 *
 * - has a fom_state_typedef
 * - has a jacobian_typedef
 * - has a jacobianCRef
 * - has a templated applyMapping(OperandType, fom_state_type)
 *
*/
template<class T, class OperandType, class ResultType = void, class Enable = void>
struct decoder : std::false_type{};

template<class T, class OperandType, class ResultType>
struct decoder<
  T, OperandType, ResultType,
  ::pressio::mpl::enable_if_t<
    mpl::not_void<ResultType>::value
    and ::pressio::has_fom_state_typedef<T>::value
    and ::pressio::has_jacobian_typedef<T>::value
    and std::is_same<ResultType, typename T::fom_state_type>::value
    and ::pressio::rom::decoder_jacobian<typename T::jacobian_type>::value
    and ::pressio::rom::has_const_get_reference_to_jacobian<T, typename T::jacobian_type>::value
    and ::pressio::rom::has_const_apply_mapping_accept_operand_result_return_void<T, OperandType, ResultType>::value
    and ::pressio::rom::has_nonconst_update_jacobian_method_accept_operand_return_void<T, OperandType>::value
    >
  > : std::true_type{};

template<class T, class OperandType>
struct decoder<
  T, OperandType, void,
  ::pressio::mpl::enable_if_t<
    ::pressio::has_fom_state_typedef<T>::value
    and ::pressio::has_jacobian_typedef<T>::value
    and ::pressio::rom::decoder_jacobian<typename T::jacobian_type>::value
    and ::pressio::rom::has_const_get_reference_to_jacobian<T, typename T::jacobian_type>::value
    and ::pressio::rom::has_const_apply_mapping_accept_operand_result_return_void<T, OperandType, typename T::fom_state_type>::value
    and ::pressio::rom::has_nonconst_update_jacobian_method_accept_operand_return_void<T, OperandType>::value
    >
  > : std::true_type{};

}}
#endif  // ROM_CONSTRAINTS_ROM_DECODER_HPP_
