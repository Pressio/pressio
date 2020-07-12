/*
//@HEADER
// ************************************************************************
//
// rom_continuous_time_system_maskable_rom.hpp
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

#ifndef rom_continuous_time_system_maskable_rom_HPP_
#define rom_continuous_time_system_maskable_rom_HPP_

namespace pressio{ namespace rom{ namespace concepts {

template<typename T, typename enable = void>
struct continuous_time_system_maskable_rom : std::false_type{};

template<typename T>
struct continuous_time_system_maskable_rom<
  T,
  mpl::enable_if_t<
    ::pressio::rom::concepts::continuous_time_system<T>::value and 
  	// createApplyMaskResult for residual   
    ::pressio::rom::predicates::has_const_create_apply_mask_result_method_accept_operand_return_result<
        T, typename T::velocity_type, typename T::velocity_type>::value and
  	// createApplyMaskResult for dense operator
    ::pressio::rom::predicates::has_const_create_apply_mask_result_method_accept_operand_return_result<
        T, typename T::dense_matrix_type, typename T::dense_matrix_type >::value and
   // applyMask for residual 
    ::pressio::rom::predicates::has_const_apply_mask_method_accept_operand_time_result_return_void<
        T, typename T::velocity_type, typename T::scalar_type, typename T::velocity_type >::value and            
  	// applyMask for dense operator
    ::pressio::rom::predicates::has_const_apply_mask_method_accept_operand_time_result_return_void<
        T, typename T::dense_matrix_type, typename T::scalar_type, typename T::dense_matrix_type >::value 
    >
  > : std::true_type{};


}}} // namespace pressio::rom::concepts
#endif
