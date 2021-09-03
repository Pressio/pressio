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

#include "../predicates/rom_has_const_apply_mapping_accept_operand_result_return_void.hpp"
#include "../predicates/rom_has_const_get_reference_to_jacobian.hpp"
#include "../predicates/rom_has_const_update_jacobian_method_accept_operand_return_void.hpp"
#include "../predicates/rom_has_nonconst_update_jacobian_method_accept_operand_return_void.hpp"
#include "../constraints/rom_decoder_jacobian.hpp"
#include "../constraints/rom_decoder.hpp"
#include "../constraints/rom_fom_types.hpp"

namespace pressio{ namespace rom{ namespace impl{

template <class FomStateType, class JacobianMatrixType>
struct LinearDecoder
{
  static_assert
  (::pressio::rom::decoder_jacobian< mpl::remove_cvref_t<JacobianMatrixType> >::value,
   "Invalid decoder's jacobian type");

  static_assert
  (::pressio::rom::fom_state< mpl::remove_cvref_t<FomStateType> >::value,
   "Invalid FomStateType");

  // these aliases are required because other ROM classes detect them
  using jacobian_type  = mpl::remove_cvref_t<JacobianMatrixType>;
  using fom_state_type = FomStateType;

public:
  LinearDecoder() = delete;
  LinearDecoder(const LinearDecoder &) = default;
  LinearDecoder & operator=(const LinearDecoder &) = default;
  LinearDecoder(LinearDecoder &&) = default;
  LinearDecoder & operator=(LinearDecoder &&) = default;
  ~LinearDecoder() = default;

  LinearDecoder(jacobian_type && matIn) 
    : jacobianOfDecoder_(std::move(matIn)){}

  LinearDecoder(const jacobian_type & matIn) 
    : jacobianOfDecoder_(matIn){}

public:
  // applyMapping must be templated because the type of the
  // generalized coordinates is not necessarily known.
  // In fact, its type is something that behaves like a vector.
  // For example, for WLS the operand here is an expression, e.g. span.

  // specialize for when gen_coords_t is a rank=1 
  template <class gen_coords_t>
  mpl::enable_if_t<::pressio::Traits<gen_coords_t>::rank == 1>
  applyMapping(const gen_coords_t & operand, fom_state_type & result) const
  {
    static_assert
      (::pressio::are_scalar_compatible<gen_coords_t, fom_state_type>::value, 
        "Types are not scalar compatible");

    using scalar_t = typename ::pressio::Traits<fom_state_type>::scalar_type;

    constexpr auto zero = ::pressio::utils::Constants<scalar_t>::zero();
    constexpr auto one  = ::pressio::utils::Constants<scalar_t>::one();
    ::pressio::ops::product(::pressio::nontranspose(), one,
			    jacobianOfDecoder_, operand, zero, result);
  }

  // specialize for when gen_coords_t is a rank=1 
  template <class gen_coords_t>
  mpl::enable_if_t<::pressio::Traits<gen_coords_t>::rank >= 2>
  applyMapping(const gen_coords_t & operand, fom_state_type & result) const
  {
    static_assert
      (::pressio::are_scalar_compatible<gen_coords_t, fom_state_type>::value, 
        "Types are not scalar compatible");

    using scalar_t = typename ::pressio::Traits<fom_state_type>::scalar_type;

    constexpr auto zero = ::pressio::utils::Constants<scalar_t>::zero();
    constexpr auto one  = ::pressio::utils::Constants<scalar_t>::one();
    ::pressio::ops::product(::pressio::nontranspose(), ::pressio::nontranspose(),
			    one, jacobianOfDecoder_, operand, zero, result);
  }

  const jacobian_type & jacobianCRef() const{
    return jacobianOfDecoder_;
  }

  template<class gen_coords_t>
  void updateJacobian(const gen_coords_t & genCoordinates)
  {
    //no op: the Jacobian is constant
  }

private:
  const jacobian_type jacobianOfDecoder_ = {};

};

}}}//end namespace pressio::rom::impl
#endif  // ROM_DECODER_IMPL_ROM_LINEAR_DECODER_PRESSIO_OPS_HPP_
