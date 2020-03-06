/*
//@HEADER
// ************************************************************************
//
// rom_is_legitimate_decoder_type.hpp
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

#ifndef ROM_IS_LEGITIMATE_DECODER_TYPE_HPP_
#define ROM_IS_LEGITIMATE_DECODER_TYPE_HPP_

namespace pressio{ namespace rom{ namespace meta {

template <typename T, typename jacobian_t, typename = void>
struct has_get_reference_to_jacobian : std::false_type{};

template <typename T, typename jacobian_t>
struct has_get_reference_to_jacobian<
  T, jacobian_t,
  mpl::enable_if_t<
    mpl::is_same<
      decltype( std::declval<T const &>().getReferenceToJacobian() ),
      const jacobian_t &
      >::value
    >
  > : std::true_type{};
//--------------------------------------------------------------


template <typename T, typename arg1_t, typename arg2_t, typename = void>
struct has_apply_mapping_two_args : std::false_type{};

template <typename T, typename arg1_t, typename arg2_t>
struct has_apply_mapping_two_args<
  T, arg1_t, arg2_t,
  mpl::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T const &>().applyMapping
       (
	std::declval<arg1_t const &>(), std::declval<arg2_t &>()
	)
       )
      >::value
    >
  > : std::true_type{};
//--------------------------------------------------------------


/*
 * A type is a legitimate decoder for LSPG if:
 *
 * - has a jacobian_typedef
 * - has a getReferenceToJacobian
 * - template applyMapping(operand_t, result_t)
 *
*/
template<
  typename T,
  typename apply_map_operand_t,
  typename apply_map_result_t,
  typename enable = void
  >
struct is_legitimate_decoder_type
  : std::false_type{};

template<
  typename T,
  typename apply_map_operand_t,
  typename apply_map_result_t
  >
struct is_legitimate_decoder_type<
  T, apply_map_operand_t, apply_map_result_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::ode::meta::has_jacobian_typedef<T>::value and
    has_get_reference_to_jacobian<T, typename T::jacobian_type>::value and
    has_apply_mapping_two_args<T, apply_map_operand_t, apply_map_result_t>::value
    >
  > : std::true_type{};

}}} // namespace pressio::rom::meta
#endif
