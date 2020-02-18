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

namespace pressio{ namespace rom{

template <
typename matrix_type, typename ops_t, 
typename rom_state_type, typename fom_state_type
>
struct PyLinearDecoder<
  matrix_type, ops_t, rom_state_type, fom_state_type,
  mpl::enable_if_t<
    ::pressio::containers::meta::is_array_pybind11<matrix_type>::value
    >
  >
  : public DecoderBase<PyLinearDecoder<matrix_type, ops_t, rom_state_type, fom_state_type>, 
  matrix_type, rom_state_type, fom_state_type>
{

  using this_t	    = PyLinearDecoder<matrix_type, ops_t, rom_state_type, fom_state_type>;
  using base_t	    = DecoderBase<this_t, matrix_type,  rom_state_type, fom_state_type>;
  using jacobian_t  = matrix_type;
  using scalar_t    = double;
  static_assert( mpl::is_same<scalar_t,
		 typename matrix_type::value_type>::value,
		 "PyLinearDecoder: Scalar types don't match");
  using rom_state_t = rom_state_type;
  using fom_state_t = fom_state_type;

private:
  friend base_t;
  matrix_type phi_ = {};

  // here we do this conditional type because it seems when ops_t= pybind11::object
  // it only works if we copy the object. Need to figure out if we can leave ptr in all cases.
  typename std::conditional<
    mpl::is_same<ops_t, pybind11::object>::value, ops_t,
    const ops_t *
    >::type ops_ = {};

  pybind11::object numpy_ = pybind11::module::import("numpy");
  pybind11::object spyblas_ = pybind11::module::import("scipy.linalg.blas");

public:
  PyLinearDecoder() = delete;

  template <typename _ops_t = ops_t, mpl::enable_if_t<std::is_void<_ops_t>::value> * = nullptr>
  PyLinearDecoder(const jacobian_t & matIn)
      // since matrix is a pybind11::array_t, use deep copy so that we own the phi here
    : phi_( jacobian_t(const_cast<jacobian_t &>(matIn).request() ))
  {}

  template <
    typename _ops_t = ops_t,
    mpl::enable_if_t< mpl::is_same<_ops_t, pybind11::object>::value > * = nullptr
    >
  PyLinearDecoder(const jacobian_t & matIn,
		  const _ops_t ops)
      // since matrix is a pybind11::array_t, use deep copy so that we own the phi here
    : phi_( jacobian_t(const_cast<jacobian_t &>(matIn).request() )),
      ops_{ops}
  {}

  ~PyLinearDecoder() = default;

private:
  const jacobian_t & getReferenceToJacobianImpl() const{
    return phi_;
  }

  /* if ops_void, and phi has col-major order, use blas*/
  template <
  typename operand_t,
  typename result_t,
  typename _ops_t = ops_t,
  typename _matrix_type = matrix_type,
  mpl::enable_if_t<
    ::pressio::containers::meta::is_fstyle_array_pybind11<_matrix_type>::value and
    std::is_void<_ops_t>::value
    > * = nullptr
  >
  void applyMappingImpl(const operand_t & operandObj,
			result_t & resultObj) const
  {
    constexpr auto zero = ::pressio::utils::constants::zero<scalar_t>();
    constexpr auto one  = ::pressio::utils::constants::one<scalar_t>();
    constexpr auto izero = ::pressio::utils::constants::zero<int>();
    constexpr auto ione  = ::pressio::utils::constants::one<int>();
    constexpr auto transA = izero;
    // overwrite y passed in to dgemv
    constexpr auto owy = ione;

    spyblas_.attr("dgemv")(one, phi_, operandObj, zero, resultObj, izero, ione, izero, ione, transA, owy);
  }

  /* if ops_t == void, phi has row-major order, use numpy*/
  template <
    typename operand_t,
    typename result_t,
    typename _ops_t = ops_t,
    typename _matrix_type = matrix_type,
    mpl::enable_if_t<
      ::pressio::containers::meta::is_cstyle_array_pybind11<_matrix_type>::value and
      std::is_void<_ops_t>::value
      > * = nullptr
  >
  void applyMappingImpl(const operand_t & operandObj,
			result_t & resultObj) const
  {
    //pybind11::object numpy = pybind11::module::import("numpy");
    // this is typically a matrix vec product. So  use matmul
    resultObj = numpy_.attr("dot")(phi_, operandObj);
  }

  /* if ops_t == pybind11::object*/
  template <
    typename operand_t,
    typename result_t,
    typename _ops_t = ops_t,
    typename _matrix_type = matrix_type,
    mpl::enable_if_t<
      ::pressio::containers::meta::is_cstyle_array_pybind11<_matrix_type>::value and
      mpl::is_same<_ops_t, pybind11::object>::value
      > * = nullptr
  >
  void applyMappingImpl(const operand_t & operandObj,
			result_t & resultObj) const
  {
    ops_.attr("multiply")(phi_, false, operandObj, false, resultObj);
  }


};//end class

}}//end namespace pressio::rom
#endif
#endif
