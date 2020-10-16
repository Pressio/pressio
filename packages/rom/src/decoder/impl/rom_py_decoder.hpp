/*
//@HEADER
// ************************************************************************
//
// rom_py_decoder.hpp
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

#ifndef ROM_DECODER_IMPL_ROM_PY_DECODER_HPP_
#define ROM_DECODER_IMPL_ROM_PY_DECODER_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <typename ...Args>
struct PyDecoder;

template <typename matrix_type, typename fom_state_type>
struct PyDecoder<
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::predicates::is_dense_matrix_wrapper_pybind11<matrix_type>::value and
    ::pressio::containers::predicates::is_vector_wrapper_pybind11<fom_state_type>::value
    >,
  matrix_type, fom_state_type
  >
{
  using jacobian_type  = matrix_type;
  enum mappingKind{ Null, Linear, Custom };

private:
  using scalar_t  =
    typename ::pressio::containers::details::traits<fom_state_type>::scalar_t;
  using jacobian_native_t =
    typename ::pressio::containers::details::traits<jacobian_type>::wrapped_t;
  using fom_native_t =
    typename ::pressio::containers::details::traits<fom_state_type>::wrapped_t;

  matrix_type mappingJacobian_ = {};
  mappingKind kind_ = {};
  pybind11::object customMapper_ = {};

public:
  PyDecoder() = delete;
  PyDecoder(const PyDecoder &) = default;
  PyDecoder & operator=(const PyDecoder &) = default;
  PyDecoder(PyDecoder &&) = default;
  PyDecoder & operator=(PyDecoder &&) = default;
  ~PyDecoder() = default;

  PyDecoder(const jacobian_native_t & jacobianMatrixIn)
    // note that we "view" the native object, we don't deep copy it.
    // if the mapping jacobian changes on the python side, it reflects here
    : mappingJacobian_(jacobianMatrixIn, ::pressio::view()),
      kind_(mappingKind::Linear)
  {}

  // here the description is not necessarily needed but it is important
  // to keep because it enables the right overload. Otherwise the interpreter
  // would pick this overload even if passing a numpy array because
  // a numpy array is also a python object.
  PyDecoder(pybind11::object customMapper, std::string description)
    // note that we "view" the native object, we don't deep copy it.
    // if the mapping jacobian changes on the python side, it reflects here
    : mappingJacobian_(customMapper.attr("jacobian")(), ::pressio::view()),
      kind_(mappingKind::Custom),
      customMapper_(customMapper)
  {}

  // applyMapping is templated because operand_t can be rom_state_type but
  // can also be an expression based on rom_state_type (e.g. for WLS)
  template <typename operand_t, typename fom_state_t = fom_state_type>
  void applyMapping(const operand_t & operand, fom_state_t & result) const
  {
    if (kind_ == mappingKind::Linear)
    {
      constexpr auto zero = ::pressio::utils::constants<scalar_t>::zero();
      constexpr auto one  = ::pressio::utils::constants<scalar_t>::one();
      ::pressio::ops::product(::pressio::nontranspose(), one,
			      mappingJacobian_, operand, zero, result);
    }
    else if(kind_ == mappingKind::Custom)
    {
      customMapper_.attr("applyMapping")(*operand.data(), *result.data());
    }
    else
      throw std::runtime_error("Invalid mapping kind enum");
  }

  const jacobian_type & jacobianCRef() const{
    return mappingJacobian_;
  }

  template<typename gen_coords_t>
  void updateJacobian(const gen_coords_t & genCoordinates) const
  {
    if (kind_ == mappingKind::Linear){
      // no op
    }
    else if(kind_ == mappingKind::Custom)
    {
      customMapper_.attr("updateJacobian")(*genCoordinates.data());
    }
    else
      throw std::runtime_error("Invalid mapping kind enum");
  }

};

}}}//end namespace pressio::rom::impl
#endif  // ROM_DECODER_IMPL_ROM_PY_DECODER_HPP_
