/*
//@HEADER
// ************************************************************************
//
// rom_linear_decoder_pressio_ops.hpp
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

#ifndef ROM_DECODER_IMPL_ROM_LINEAR_DECODER_PRESSIO_OPS_HPP_
#define ROM_DECODER_IMPL_ROM_LINEAR_DECODER_PRESSIO_OPS_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <typename jacobian_matrix_type, typename fom_state_type>
struct LinearDecoderWithPressioOps
{
  using jacobian_type  = jacobian_matrix_type;
  using scalar_t =
    typename ::pressio::containers::details::traits<fom_state_type>::scalar_t;
  using jacobian_native_t =
    typename ::pressio::containers::details::traits<jacobian_type>::wrapped_t;

private:
  const jacobian_type jacobianOfDecoder_ = {};

public:
  LinearDecoderWithPressioOps() = delete;

  // /* note that here we make the constructor templated
  //    such that we can pass either
  //    the pressio wrapper type 'jacobian_matrix_type' or
  //    the native type wrapped by it 'jacobian_native_t'.
  //    Also, by templating this constructor we enable
  //    universal reference so that it is forwarded accordingly.
  //  */
  template<typename T>
  LinearDecoderWithPressioOps(T && matIn)
    : jacobianOfDecoder_(std::forward<T>(matIn)){}

public:
  // applyMapping is templated because the type of the generalized coordinates
  // is something that behaves like a vector but not strictly a vector:
  // for LSPG it is a "concrete" pressio::Vector but for WLS is an expression
  template <typename gen_coords_t, typename fom_state_t = fom_state_type>
  void applyMapping(const gen_coords_t & operand, fom_state_t & result) const
  {
    constexpr auto zero = ::pressio::utils::constants<scalar_t>::zero();
    constexpr auto one  = ::pressio::utils::constants<scalar_t>::one();
    ::pressio::ops::product(::pressio::nontranspose(), one,
			    jacobianOfDecoder_, operand, zero, result);
  }

  const jacobian_type & getReferenceToJacobian() const{
    return jacobianOfDecoder_;
  }

  // template<typename gen_coords_t>
  // void updateJacobian(const gen_coords_t & genCoordinates) const
  // {
  //   //no op
  // }
};

}}}//end namespace pressio::rom::impl
#endif  // ROM_DECODER_IMPL_ROM_LINEAR_DECODER_PRESSIO_OPS_HPP_
