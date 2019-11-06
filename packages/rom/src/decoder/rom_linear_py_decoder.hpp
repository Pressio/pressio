/*
//@HEADER
// ************************************************************************
//
// rom_linear_py_decoder.hpp
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

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
#ifndef ROM_LINEAR_PY_DECODER_HPP_
#define ROM_LINEAR_PY_DECODER_HPP_

#include "rom_decoder_base.hpp"
#include "../rom_fwd.hpp"

namespace pressio{ namespace rom{

template <
  typename matrix_type,
  typename ops_t
  >
struct PyLinearDecoder<
  matrix_type, ops_t,
  mpl::enable_if_t<
    ::pressio::containers::meta::is_array_pybind11<matrix_type>::value
    >
  >
  : public DecoderBase<
  PyLinearDecoder<matrix_type, ops_t>, matrix_type
  >
{

  using this_t	    = PyLinearDecoder<matrix_type, ops_t>;
  using base_t	    = DecoderBase<this_t, matrix_type>;
  using jacobian_t  = matrix_type;
  using scalar_t    = double;
  static_assert( mpl::is_same<scalar_t,
		 typename matrix_type::value_type>::value,
		 "PyLinearDecoder: Scalar types don't match");

private:
  friend base_t;
  matrix_type phi_ = {};
  ops_t ops_ = {};

public:
  PyLinearDecoder() = delete;

  // since matrix is a pybind11::array_t, the following
  // uses view semantics, so NOT a deep copy. Hence, the object
  // is owened on the python side. This is fine, beucase this is the
  // only object who owns the basis vectors. if we want to do a deep copy,
  // then we have to init phi_ using the buffer of matIn like we do
  // in ode_storage for example
  PyLinearDecoder(const jacobian_t & matIn,
		  const ops_t ops)
    : phi_( jacobian_t(const_cast<jacobian_t &>(matIn).request() )),
      ops_{ops}
  {}

  ~PyLinearDecoder() = default;

  template <typename operand_t, typename result_t>
  void _applyMappingTest(const operand_t & operandObj,
			result_t & resultObj) const{
    ops_.attr("multiply")(phi_, false, operandObj, false, resultObj);
  }

protected:
  template <typename operand_t, typename result_t>
  void applyMappingImpl(const operand_t & operandObj,
			result_t & resultObj) const{
    ops_.attr("multiply")(phi_, false, operandObj, false, resultObj);
  }

  const jacobian_t & getReferenceToJacobianImpl() const{
    return phi_;
  }

};//end class

}}//end namespace pressio::rom
#endif
#endif
