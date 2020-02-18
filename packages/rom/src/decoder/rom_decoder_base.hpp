/*
//@HEADER
// ************************************************************************
//
// rom_decoder_base.hpp
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

#ifndef ROM_DECODER_BASE_HPP_
#define ROM_DECODER_BASE_HPP_

namespace pressio{ namespace rom{

template <
  typename derived_type,
  typename jac_matrix_type,
  typename rom_state_type,
  typename fom_state_type
  >
struct DecoderBase
{
  using this_t = DecoderBase<derived_type,
			     jac_matrix_type,
			     rom_state_type,
			     fom_state_type>;

  // template <typename operand_t, typename result_t>
  // void applyMapping(const operand_t & operandObj,
  // 		    result_t & result) const  {
  //   static_cast<const derived_type &>(*this).applyMappingImpl(operandObj, result);
  // }

  template <typename operand_t>
  void applyMapping(const operand_t & operandObj, fom_state_type & result) const
  {
    static_cast<const derived_type &>(*this).template applyMappingImpl<operand_t>(operandObj, result);
  }

  const jac_matrix_type & getReferenceToJacobian() const {
    return static_cast<const derived_type &>(*this).getReferenceToJacobianImpl();
  }

  DecoderBase() = default;
  ~DecoderBase() = default;

};//end class

}}//end namespace pressio::rom
#endif
