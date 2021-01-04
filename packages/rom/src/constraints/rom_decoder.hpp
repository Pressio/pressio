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

namespace pressio{ namespace rom{ namespace constraints {

/*
 * A type T is a legitimate decoder if:
 *
 * - has a fom_state_typedef
 * - has a jacobian_typedef
 * - has a jacobianCRef
 * - template applyMapping(operand_t, fom_state_type)
 *
*/
template<
  typename T,
  typename operand_t,
  typename result_t = void,
  typename enable = void
  >
struct decoder : std::false_type{};

template<
  typename T,
  typename operand_t,
  typename result_t
  >
struct decoder<
  T, operand_t, result_t,
  ::pressio::mpl::enable_if_t<
    mpl::not_void<result_t>::value
    and ::pressio::rom::predicates::has_fom_state_typedef<T>::value
    and ::pressio::ode::predicates::has_jacobian_typedef<T>::value
    and std::is_same<result_t, typename T::fom_state_type>::value
    and ::pressio::rom::constraints::decoder_jacobian<typename T::jacobian_type>::value
    and
    ::pressio::rom::predicates::has_const_get_reference_to_jacobian<
      T, typename T::jacobian_type>::value
    and
    ::pressio::rom::predicates::has_const_apply_mapping_accept_operand_result_return_void<
      T, operand_t, result_t>::value
    and
    ::pressio::rom::predicates::has_nonconst_update_jacobian_method_accept_operand_return_void<
      T, operand_t>::value
    >
  > : std::true_type{};

template<typename T, typename operand_t>
struct decoder<
  T, operand_t, void,
  ::pressio::mpl::enable_if_t<
    ::pressio::rom::predicates::has_fom_state_typedef<T>::value
    and ::pressio::ode::predicates::has_jacobian_typedef<T>::value
    and ::pressio::rom::constraints::decoder_jacobian<typename T::jacobian_type>::value
    and
    ::pressio::rom::predicates::has_const_get_reference_to_jacobian<
      T, typename T::jacobian_type>::value
    and
    ::pressio::rom::predicates::has_const_apply_mapping_accept_operand_result_return_void<
      T, operand_t, typename T::fom_state_type>::value
    and
    ::pressio::rom::predicates::has_nonconst_update_jacobian_method_accept_operand_return_void<
      T, operand_t>::value
    >
  > : std::true_type{};

}}} // namespace pressio::rom::constraints
#endif  // ROM_CONSTRAINTS_ROM_DECODER_HPP_
