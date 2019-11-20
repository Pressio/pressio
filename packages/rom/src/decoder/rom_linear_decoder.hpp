/*
//@HEADER
// ************************************************************************
//
// rom_linear_decoder.hpp
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

#ifndef ROM_LINEAR_DECODER_HPP_
#define ROM_LINEAR_DECODER_HPP_

#include "rom_decoder_base.hpp"
#include "../rom_fwd.hpp"

namespace pressio{ namespace rom{

template <
  typename matrix_type,
  template <typename...> class
  wrapper_operator_t = MultiVectorOperator
  >
struct LinearDecoder
  : public DecoderBase<
  LinearDecoder<matrix_type, wrapper_operator_t>,
  matrix_type>{

  using this_t	    = LinearDecoder<matrix_type, wrapper_operator_t>;
  using base_t	    = DecoderBase<this_t, matrix_type>;
  using matrix_op_t = wrapper_operator_t<matrix_type>;
  using jacobian_t  = matrix_type;

private:
  friend base_t;
  matrix_op_t phi_ = {};

public:
  LinearDecoder() = delete;

  LinearDecoder(const jacobian_t & matIn)
    : phi_(matIn){}

  ~LinearDecoder() = default;

private:
  template <typename operand_t, typename result_t>
  void applyMappingImpl(const operand_t & operandObj,
			result_t & resultObj) const{
    phi_.apply(operandObj, resultObj);
  }

  const jacobian_t & getReferenceToJacobianImpl() const{
    return phi_.getRefToOperator();
  }

};//end class

}}//end namespace pressio::rom
#endif
