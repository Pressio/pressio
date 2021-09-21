/*
//@HEADER
// ************************************************************************
//
// rom_masker.hpp
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

#ifndef ROM_CONSTRAINTS_ROM_MASKER_EXPLICIT_HPP_
#define ROM_CONSTRAINTS_ROM_MASKER_EXPLICIT_HPP_

namespace pressio{ namespace rom{

template<class T, class TimeType, class FomVelocityType, class enable = void>
struct masker_explicit_galerkin : std::false_type{};

template<class T, class TimeType, class FomVelocityType>
struct masker_explicit_galerkin<
  T, TimeType, FomVelocityType,
  mpl::enable_if_t<
    // createApplyMaskResult for FomVelocityType
    ::pressio::rom::has_const_create_apply_mask_result_method_accept_operand_return_result<
      T, FomVelocityType, FomVelocityType>::value
    and
    // callable for FomVelocityType
    std::is_same<
      decltype(std::declval<T const>()
	       (
		std::declval<FomVelocityType const&>(),
		std::declval<TimeType>(),
		std::declval<FomVelocityType &>()
		)
	       ), void
      >::value
    >
  > : std::true_type{};

template<
  class T,
  class TimeType,
  class Operand1Type,
  class Operand2Type,
  class enable = void
  >
struct masker_implicit_galerkin : std::false_type{};

template<
  class T,
  class TimeType,
  class Operand1Type,
  class Operand2Type
  >
struct masker_implicit_galerkin<
  T, TimeType, Operand1Type, Operand2Type,
  mpl::enable_if_t<
    // createApplyMaskResult for Operand1Type
    ::pressio::rom::has_const_create_apply_mask_result_method_accept_operand_return_result<
      T, Operand1Type, Operand1Type>::value
    and
    // createApplyMaskResult for Operand2Type
    ::pressio::rom::has_const_create_apply_mask_result_method_accept_operand_return_result<
      T, Operand2Type, Operand2Type>::value
    and
    // callable for operand1
    std::is_same<
      decltype(std::declval<T const>()
	       (
		std::declval<Operand1Type const&>(),
		std::declval<TimeType>(),
		std::declval<Operand1Type &>()
		)
	       ), void
      >::value
    and
    // callable for operand1
    std::is_same<
      decltype(std::declval<T const>()
	       (
		std::declval<Operand2Type const&>(),
		std::declval<TimeType>(),
		std::declval<Operand2Type &>()
		)
	       ), void
      >::value
   >
  > : std::true_type{};

}}
#endif  // ROM_GALERKIN_CONSTRAINTS_ROM_MASKER_EXPLICIT_HPP_
