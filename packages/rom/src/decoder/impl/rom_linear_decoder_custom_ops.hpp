/*
//@HEADER
// ************************************************************************
//
// rom_linear_decoder_custom_ops.hpp
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

#ifndef ROM_DECODER_IMPL_ROM_LINEAR_DECODER_CUSTOM_OPS_HPP_
#define ROM_DECODER_IMPL_ROM_LINEAR_DECODER_CUSTOM_OPS_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <typename matrix_type, typename fom_state_type, typename ops_t>
struct LinearDecoderWithCustomOps
{
  using jacobian_type  = matrix_type;
  using scalar_t =
    typename ::pressio::containers::details::traits<fom_state_type>::scalar_t;
  using jacobian_native_t =
    typename ::pressio::containers::details::traits<matrix_type>::wrapped_t;

private:
  const matrix_type mappingJacobian_ = {};
  std::reference_wrapper<const ops_t> udOps_;

public:
  LinearDecoderWithCustomOps() = delete;
  LinearDecoderWithCustomOps(const LinearDecoderWithCustomOps &) = default;
  LinearDecoderWithCustomOps & operator=(const LinearDecoderWithCustomOps &) = default;
  LinearDecoderWithCustomOps(LinearDecoderWithCustomOps &&) = default;
  LinearDecoderWithCustomOps & operator=(LinearDecoderWithCustomOps &&) = default;
  ~LinearDecoderWithCustomOps() = default;


  // /* note that here we make the constructor templated
  //    such that we can pass either
  //    the pressio wrapper type 'jacobian_matrix_type' or
  //    the native type wrapped by it 'jacobian_native_t'.
  //    Also, by templating this constructor we enable
  //    universal reference so that it is forwarded accordingly.
  //  */
  template<typename T>
  LinearDecoderWithCustomOps(T && matIn,
			     const ops_t & udOps)
    : mappingJacobian_(std::forward<T>(matIn)),
      udOps_{udOps}{}

  // applyMapping is templated because the type of the generalized coordinates
  // is something that behaves like a vector but not strictly a vector:
  // for LSPG it is a "concrete" pressio::Vector but for WLS is an expression
  template <typename gen_coords_t>
  void applyMapping(const gen_coords_t & operand, fom_state_type & result) const
  {
    constexpr auto zero = ::pressio::utils::constants<scalar_t>::zero();
    constexpr auto one  = ::pressio::utils::constants<scalar_t>::one();
    udOps_.get().product(::pressio::nontranspose(), one, *mappingJacobian_.data(),
		   operand, zero, *result.data());
  }

  const jacobian_type & jacobianCRef() const{
    return mappingJacobian_;
  }

  template<typename gen_coords_t>
  void updateJacobian(const gen_coords_t & genCoordinates) const
  {
    //no op
  }

};//end

}}}//end namespace pressio::rom::impl
#endif  // ROM_DECODER_IMPL_ROM_LINEAR_DECODER_CUSTOM_OPS_HPP_
