/*
//@HEADER
// ************************************************************************
//
// rom_create_collocation_projector.hpp
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

#ifndef ROM_GALERKIN_ROM_CREATE_COLLOCATION_PROJECTOR_HPP_
#define ROM_GALERKIN_ROM_CREATE_COLLOCATION_PROJECTOR_HPP_

namespace pressio{ namespace rom{ namespace galerkin{

template<typename decoder_type, typename collocator_t>
impl::ArbitraryProjector<typename decoder_type::jacobian_type, void>
createCollocationProjector(const decoder_type & decoder,
			   collocator_t && collocator)
{
  // call the collocator passing the native jacobian object
  const auto & nativeDecoderJacobian = *decoder.jacobianCRef().data();
  auto matrix = collocator.sampleRows(nativeDecoderJacobian);

  // this is the wrapped type for the decoder's jacobian
  using dec_jac_t = typename decoder_type::jacobian_type;

  // to construct the arbitrary projector, move the matrix
  // so that we hopefully don't waste a copy
  using return_t = impl::ArbitraryProjector<dec_jac_t, void>;
  return return_t(std::move(matrix));
}

template<typename decoder_type, typename collocator_t, typename ops_t>
impl::ArbitraryProjector<typename decoder_type::jacobian_type, ops_t>
createCollocationProjector(const decoder_type & decoder,
			   collocator_t && collocator,
			   const ops_t & udOps)
{
  // call the collocator passing the native jacobian object
  const auto & nativeDecoderJacobian = *decoder.jacobianCRef().data();
  auto matrix = collocator.sampleRows(nativeDecoderJacobian);

  // this is the wrapped type for the decoder's jacobian
  using dec_jac_t = typename decoder_type::jacobian_type;

  // to construct the arbitrary projector, move the matrix
  // so that we hopefully don't waste a copy
  using return_t = impl::ArbitraryProjector<dec_jac_t, ops_t>;
  return return_t(std::move(matrix), udOps);
}

}}}//end namespace pressio::rom::galerkin
#endif  // ROM_GALERKIN_ROM_CREATE_COLLOCATION_PROJECTOR_HPP_
