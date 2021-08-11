/*
//@HEADER
// ************************************************************************
//
// rom_linear_decoder_specializer.hpp
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

#ifndef ROM_DECODER_IMPL_ROM_LINEAR_DECODER_SPECIALIZER_HPP_
#define ROM_DECODER_IMPL_ROM_LINEAR_DECODER_SPECIALIZER_HPP_

#include "rom_linear_decoder_pressio_ops.hpp"
#include "rom_linear_decoder_custom_ops.hpp"

namespace pressio{ namespace rom{ namespace impl{

template <typename ... Args>
struct LinearDecoderSpecializer;

// matrix_type, fom_state_type
template <typename matrix_type, typename fom_state_type>
struct LinearDecoderSpecializer<
  void, matrix_type, fom_state_type
  >
{
  static_assert
  (::pressio::rom::constraints::decoder_jacobian<matrix_type>::value,
   "Template arg passed to the LinearDecoder not a valid decoder's Jacobian type.");

  static_assert
  (::pressio::containers::predicates::is_wrapper<fom_state_type>::value,
   "Fom state type template arg passed to the LinearDecoder class must be a wrapper");

  using type = LinearDecoderWithPressioOps<matrix_type, fom_state_type>;
};

// matrix_type, fom_state_type, ud_ops
template <typename matrix_type, typename fom_state_type, typename ud_ops_t>
struct LinearDecoderSpecializer<
  void, matrix_type, fom_state_type, ud_ops_t
  >
{
  static_assert
  (::pressio::rom::constraints::decoder_jacobian<matrix_type>::value,
   "Template arg passed to the LinearDecoder not a valid decoder's Jacobian type.");

  static_assert
  (::pressio::containers::predicates::is_wrapper<fom_state_type>::value,
   "Fom state type template arg passed to the LinearDecoder class must be a wrapper");

  using type = LinearDecoderWithCustomOps<matrix_type, fom_state_type,ud_ops_t>;
};

}}}//end namespace pressio::rom::impl
#endif  // ROM_DECODER_IMPL_ROM_LINEAR_DECODER_SPECIALIZER_HPP_
