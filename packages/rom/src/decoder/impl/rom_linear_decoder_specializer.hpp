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

template <
  typename matrix_type, typename rom_state_type, typename fom_state_type,
  typename ...>
struct LinearDecoderSpecializer;


template <
  typename matrix_type, typename rom_state_type, typename fom_state_type
  >
struct LinearDecoderSpecializer<
  matrix_type, rom_state_type, fom_state_type
  >
{
  /* ------------------------------------------------ */
  // the matrix type passed must be either a multivector wrapper
  // or a matrix wrapper, so we need to detect if either type if found
  static constexpr auto vMV =
    ::pressio::containers::predicates::is_multi_vector_wrapper<matrix_type>::value;
  static constexpr auto vMat =
    ::pressio::containers::predicates::is_matrix_wrapper<matrix_type>::value;
  static_assert(vMV or vMat,
		"The matrix template arg passed to the LinearDecoder class \
must be a multi-vector or matrix wrapper type. ");

  /* ------------------------------------------------ */
  // check for a valid ROM state type
  static_assert(::pressio::containers::predicates::is_vector_wrapper<rom_state_type>::value,
		"The rom state type template arg passed to the LinearDecoder class \
must be (for now) a pressio vector wrapper.");

  /* ------------------------------------------------ */
  // check for a valid FOM state type
  static_assert(::pressio::containers::predicates::is_vector_wrapper<fom_state_type>::value,
		"The fom state type template arg passed to the LinearDecoder class \
must be (for now) a pressio vector wrapper.");

  /* ------------------------------------------------ */
  using type = ::pressio::rom::impl::LinearDecoderWithPressioOps<matrix_type, rom_state_type, fom_state_type>;
};


template <
  typename matrix_type,
  typename rom_state_type,
  typename fom_state_type,
  typename ops_t
  >
struct LinearDecoderSpecializer<
  matrix_type, rom_state_type, fom_state_type, ops_t
  >
{

  /* ------------------------------------------------ */
  // the matrix type passed must be either a multivector wrapper
  // or a matrix wrapper, so we need to detect if either type if found
  static constexpr auto vMV =
    ::pressio::containers::predicates::is_multi_vector_wrapper<matrix_type>::value;
  static constexpr auto vMat =
    ::pressio::containers::predicates::is_matrix_wrapper<matrix_type>::value;
  static_assert(vMV or vMat,
		"The matrix template arg passed to the LinearDecoder class \
must be a multi-vector or matrix wrapper type. ");

  /* ------------------------------------------------ */
  // check for a valid ROM state type
  static_assert(::pressio::containers::predicates::is_vector_wrapper<rom_state_type>::value,
		"The rom state type template arg passed to the LinearDecoder class \
must be (for now) a pressio vector wrapper.");

  /* ------------------------------------------------ */
  // check for a valid FOM state type
  static_assert(::pressio::containers::predicates::is_vector_wrapper<fom_state_type>::value,
		"The fom state type template arg passed to the LinearDecoder class \
must be (for now) a pressio vector wrapper.");

  /* ------------------------------------------------ */
  // detect for valid ops
  static_assert(::pressio::rom::concepts::custom_ops_for_linear_decoder<
		ops_t, matrix_type, rom_state_type, fom_state_type
		>::value,
		"You are tring to create linear decoder with custom ops. \
The template arg passed to the LinearDecoder class representing the custom ops \
does not have an admissible API for the linear decoder class.");

  /* ------------------------------------------------ */
  using type = ::pressio::rom::impl::LinearDecoderWithCustomOps<matrix_type,
								rom_state_type,
								fom_state_type,
								ops_t>;
};

}}}//end namespace pressio::rom::impl
#endif  // ROM_DECODER_IMPL_ROM_LINEAR_DECODER_SPECIALIZER_HPP_
