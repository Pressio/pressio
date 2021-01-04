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

template <typename jacobian_matrix_type, typename fom_state_t>
struct LinearDecoderWithPressioOps
{
  static_assert
  (::pressio::rom::constraints::decoder_jacobian<jacobian_matrix_type>::value,
   "Invalid decoder's jacobian type");
  static_assert
  (::pressio::containers::predicates::is_wrapper<fom_state_t>::value,
   "Invalid fom_state type passed to decoder");

  // these aliases must be here because ROM classes detect them
  using jacobian_type  = jacobian_matrix_type;
  using fom_state_type = fom_state_t;

private:
  const jacobian_type jacobianOfDecoder_ = {};

public:
  LinearDecoderWithPressioOps() = delete;
  LinearDecoderWithPressioOps(const LinearDecoderWithPressioOps &) = default;
  LinearDecoderWithPressioOps & operator=(const LinearDecoderWithPressioOps &) = default;
  LinearDecoderWithPressioOps(LinearDecoderWithPressioOps &&) = default;
  LinearDecoderWithPressioOps & operator=(LinearDecoderWithPressioOps &&) = default;
  ~LinearDecoderWithPressioOps() = default;

  /* the constructor is templated such that we can pass either
     the pressio wrapper type 'jacobian_matrix_type' or
     the native type wrapped by it 'jacobian_native_t'.
     By templating this constructor we enable universal references
     so that it is forwarded accordingly.
   */
  template<typename T>
  LinearDecoderWithPressioOps(T && matIn)
    : jacobianOfDecoder_(std::forward<T>(matIn))
  {
    // PRESSIOLOG_DEBUG
    //   ("cnstr: possibly allocating matrix with size = ({},{}), addr = {}",
    //    &jacobianOfDecoder_,
    //    jacobianOfDecoder_.extent(0), jacobianOfDecoder_.extent(1));
  }

public:
  // applyMapping must be templated because the type of the
  // generalized coordinates is not necessarily know.
  // In fact, its type is something that behaves like a vector
  // but not strictly a vector type.
  // For example, for WLS the operand here is a pressio expression, e.g. span.
  template <typename gen_coords_t>
  mpl::enable_if_t<::pressio::containers::details::traits<gen_coords_t>::rank == 1>
  applyMapping(const gen_coords_t & operand, fom_state_type & result) const
  {
    static_assert
      (::pressio::containers::predicates::are_scalar_compatible<
       gen_coords_t, fom_state_type>::value, "Types are not scalar compatible");
    using scalar_t = typename ::pressio::containers::details::traits<
      fom_state_type>::scalar_t;

    constexpr auto zero = ::pressio::utils::constants<scalar_t>::zero();
    constexpr auto one  = ::pressio::utils::constants<scalar_t>::one();
    ::pressio::ops::product(::pressio::nontranspose(), one,
			    jacobianOfDecoder_, operand, zero, result);
  }

  template <typename gen_coords_t>
  mpl::enable_if_t<::pressio::containers::details::traits<gen_coords_t>::rank >= 2>
  applyMapping(const gen_coords_t & operand, fom_state_type & result) const
  {
    static_assert
      (::pressio::containers::predicates::are_scalar_compatible<
       gen_coords_t, fom_state_type>::value, "Types are not scalar compatible");
    using scalar_t = typename ::pressio::containers::details::traits<
      fom_state_type>::scalar_t;

    constexpr auto zero = ::pressio::utils::constants<scalar_t>::zero();
    constexpr auto one  = ::pressio::utils::constants<scalar_t>::one();
    ::pressio::ops::product(::pressio::nontranspose(), ::pressio::nontranspose(),
			    one, jacobianOfDecoder_, operand, zero, result);
  }

  const jacobian_type & jacobianCRef() const{
    return jacobianOfDecoder_;
  }

  template<typename gen_coords_t>
  void updateJacobian(const gen_coords_t & genCoordinates)
  {
    //no op: the Jacobian is constant
  }
};

}}}//end namespace pressio::rom::impl
#endif  // ROM_DECODER_IMPL_ROM_LINEAR_DECODER_PRESSIO_OPS_HPP_
