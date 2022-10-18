/*
//@HEADER
// ************************************************************************
//
// ops_elementwise_multiply.hpp
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

#ifndef OPS_TPETRA_BLOCK_OPS_ELEMENTWISE_MULTIPLY_HPP_
#define OPS_TPETRA_BLOCK_OPS_ELEMENTWISE_MULTIPLY_HPP_

namespace pressio{ namespace ops{

//----------------------------------------------------------------------
// computing elementwise:  y = beta * y + alpha * x * z
//----------------------------------------------------------------------
template <typename T, typename T1, typename T2, class alpha_t, class beta_t>
::pressio::mpl::enable_if_t<
  ::pressio::is_vector_tpetra_block<T>::value and
  ::pressio::is_vector_tpetra_block<T1>::value and
  ::pressio::is_vector_tpetra_block<T2>::value
  >
elementwise_multiply(const alpha_t & alpha,
		     const T & x,
		     const T1 & z,
		     const beta_t & beta,
		     T2 & y)
{
  assert(extent(x,0)==extent(z,0));
  assert(extent(z,0)==extent(y,0));

  auto x_tpetraview = const_cast<T &>(x).getVectorView();
  auto z_tpetraview = const_cast<T1 &>(z).getVectorView();
  auto y_tpetraview = y.getVectorView();
  y_tpetraview.elementWiseMultiply(alpha, x_tpetraview, z_tpetraview, beta);
}

}}//end namespace pressio::ops
#endif  // OPS_TPETRA_BLOCK_OPS_ELEMENTWISE_MULTIPLY_HPP_
