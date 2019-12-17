/*
//@HEADER
// ************************************************************************
//
// containers_mvec_dot_mvec.hpp
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
#ifndef CONTAINERS_SRC_OPS_EIGEN_MULTI_VECTOR_DOT_MVEC_HPP_
#define CONTAINERS_SRC_OPS_EIGEN_MULTI_VECTOR_DOT_MVEC_HPP_

#include "../containers_ops_meta.hpp"
#include "../../multi_vector/containers_multi_vector_meta.hpp"

namespace pressio{ namespace containers{ namespace ops{

/*
 * multi_vector dot multi_vector
 */

/* ----------------------------------------------
 * result_t = dense dynamic eigen matrix wrapper
 ---------------------------------------------- */
template <
  typename mvec_t,
  typename result_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_eigen<mvec_t>::value and
    ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<result_t>::value and
    ::pressio::containers::details::traits<result_t>::is_dynamic and
    ::pressio::containers::meta::wrapper_pair_have_same_scalar<mvec_t, result_t>::value
    > * = nullptr
  >
void dot(const mvec_t & A, const mvec_t & B, result_t & C)
{
  const auto nAcols = A.numVectors();
  const auto nArows = A.length();
  const auto nBcols = B.numVectors();
  const auto nBrows = B.length();
  assert(nArows == nBrows);

  // since C is dynamic, I can resize if needed
  if(C.rows() != nAcols || C.cols() != nBcols)
    C.data()->resize( nAcols, nBcols );

  *C.data() = A.data()->transpose() * (*B.data());
}

template <
  typename mvec_t,
  typename result_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_eigen<mvec_t>::value and
    ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<result_t>::value and
    ::pressio::containers::details::traits<result_t>::is_dynamic
    > * = nullptr
  >
result_t dot(const mvec_t & A, const mvec_t & B)
{
  static_assert( ::pressio::containers::meta::wrapper_pair_have_same_scalar<mvec_t, result_t>::value,
		 "MV dot MV for Eigen wrappers: the MV and result type need to have matching scalar types");

  // this is A^T B
  const auto numVecsA = A.numVectors();
  const auto numVecsB = B.numVectors();
  result_t C(numVecsA, numVecsB);
  dot(A, B, C);
  return C;
}



/* ----------------------------------------------
 * result_t = an expression
 ---------------------------------------------- */
template <
  typename mvec_t,
  typename expr_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_eigen<mvec_t>::value and
    ::pressio::containers::meta::is_expression<expr_t>::value
    > * = nullptr
  >
void dot(const mvec_t & A, const mvec_t & B, expr_t & exprObj)
{
  using mv_scalar_t   = typename ::pressio::containers::details::traits<mvec_t>::scalar_t;
  using expr_scalar_t = typename ::pressio::containers::details::traits<expr_t>::scalar_t;
  static_assert( std::is_same<mv_scalar_t, expr_scalar_t>::value,
		 "MV dot MV for Eigen wrappers: the MV and expr types need to have matching scalar types");

  const auto nAcols = A.numVectors();
  const auto nArows = A.length();
  const auto nBcols = B.numVectors();
  const auto nBrows = B.length();
  assert(nArows == nBrows);

  assert(exprObj.rows()==nAcols and exprObj.cols()==nBcols);
  exprObj() = A.data()->transpose() * (*B.data());
}

}}}//end namespace pressio::containers::ops
#endif
