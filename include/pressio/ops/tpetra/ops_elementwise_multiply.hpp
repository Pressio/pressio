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

#ifndef OPS_TPETRA_OPS_ELEMENTWISE_MULTIPLY_HPP_
#define OPS_TPETRA_OPS_ELEMENTWISE_MULTIPLY_HPP_

namespace pressio{ namespace ops{

//----------------------------------------------------------------------
// computing elementwise:  y = beta * y + alpha * x * z
//----------------------------------------------------------------------
template <typename T, typename T1, typename T2, class alpha_t, class beta_t>
::pressio::mpl::enable_if_t<
  // TPL/container specific
     (   ::pressio::is_vector_tpetra<T>::value
      || ::pressio::is_expression_column_acting_on_tpetra<T>::value)
  && (   ::pressio::is_vector_tpetra<T1>::value
      || ::pressio::is_expression_column_acting_on_tpetra<T1>::value)
  && (   ::pressio::is_vector_tpetra<T2>::value
      || ::pressio::is_expression_column_acting_on_tpetra<T2>::value)
  && Traits<T>::rank == 1
  && Traits<T1>::rank == 1
  && Traits<T2>::rank == 1
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<T, T1, T2>::value
  && std::is_convertible<alpha_t, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<beta_t,  typename ::pressio::Traits<T>::scalar_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<T>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<T>::scalar_type>::value)
  >
elementwise_multiply(const alpha_t & alpha,
		     const T & xin,
		     const T1 & zin,
		     const beta_t & beta,
		     T2 & yin)
{
  auto x = impl::get_native(xin);
  auto z = impl::get_native(zin);
  auto y = impl::get_native(yin);

  assert(::pressio::ops::extent(x, 0)==::pressio::ops::extent(z, 0));
  assert(::pressio::ops::extent(z, 0)==::pressio::ops::extent(y, 0));

  using sc_t = typename ::pressio::Traits<T>::scalar_type;
  const sc_t alpha_(alpha);
  const sc_t beta_(beta);
  y.elementWiseMultiply(alpha_, x, z, beta_);
}

}}//end namespace pressio::ops
#endif  // OPS_TPETRA_OPS_ELEMENTWISE_MULTIPLY_HPP_
