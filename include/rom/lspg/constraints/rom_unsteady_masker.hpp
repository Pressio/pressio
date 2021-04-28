/*
//@HEADER
// ************************************************************************
//
// rom_unsteady_masker.hpp
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

#ifndef ROM_LSPG_CONSTRAINTS_ROM_UNSTEADY_MASKER_HPP_
#define ROM_LSPG_CONSTRAINTS_ROM_UNSTEADY_MASKER_HPP_

namespace pressio{ namespace rom{ namespace lspg{namespace constraints {

template<
  typename T,
  typename scalar_t,
  typename operand1_t,
  typename operand2_t,
  typename enable = void
  >
struct unsteady_masker : std::false_type{};

template<
  typename T,
  typename scalar_t,
  typename operand1_t,
  typename operand2_t
  >
struct unsteady_masker<
  T, scalar_t, operand1_t, operand2_t,
  mpl::enable_if_t<
    // createApplyMaskResult for operand1_t
    ::pressio::rom::predicates::has_const_create_apply_mask_result_method_accept_operand_return_result<
      T, operand1_t, operand1_t>::value and
    // createApplyMaskResult for operand2_t
    ::pressio::rom::predicates::has_const_create_apply_mask_result_method_accept_operand_return_result<
      T, operand2_t, operand2_t>::value and
    // applyMask for operand1_t
    ::pressio::rom::predicates::has_const_apply_mask_method_accept_operand_time_result_return_void<
      T, operand1_t, scalar_t, operand1_t >::value and
    // applyMask for operand2_t
    ::pressio::rom::predicates::has_const_apply_mask_method_accept_operand_time_result_return_void<
      T, operand2_t, scalar_t, operand2_t >::value
    >
  > : std::true_type{};


// template <typename T>
// struct find_discrepancies_with_continuous_time_system_with_user_provided_apply_jacobian_maskable_api
// {
//   static_assert
//   (::pressio::rom::find_discrepancies_with_continuous_time_system_with_user_provided_apply_jacobian_api<T>::value,"");

//   static_assert
//     (::pressio::rom::predicates::has_const_create_apply_mask_result_method_accept_operand_return_result<
//      T, typename T::velocity_type, typename T::velocity_type >::value,
//      "Your continuous-time adapter class is without (or has a wrong) create apply mask to velocity method");

//   static_assert
//     (::pressio::rom::predicates::has_const_create_apply_mask_result_method_accept_operand_return_result<
//      T, typename T::dense_matrix_type, typename T::dense_matrix_type >::value,
//      "Your continuous-time adapter class is without (or has a wrong) create apply mask to dense matrix method");

//   static_assert
//     (::pressio::rom::predicates::has_const_apply_mask_method_accept_operand_time_result_return_void<
//      T, typename T::state_type, typename T::scalar_type, typename T::velocity_type >::value,
//      "Your continuous-time adapter class is without (or has a wrong) apply mask to velocity method");

//   static_assert
//     (::pressio::rom::predicates::has_const_apply_mask_method_accept_operand_time_result_return_void<
//      T, typename T::dense_matrix_type, typename T::scalar_type, typename T::dense_matrix_type >::value,
//      "Your continuous-time adapter class is without (or has a wrong) apply mask to dense matrix method");

//   static constexpr bool value = true;
// };

}}}} 
#endif  // ROM_LSPG_CONSTRAINTS_ROM_UNSTEADY_MASKER_HPP_
