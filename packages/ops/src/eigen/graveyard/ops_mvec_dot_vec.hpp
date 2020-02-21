/*
//@HEADER
// ************************************************************************
//
// containers_mvec_dot_vec.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
]//
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
#ifndef CONTAINERS_SRC_OPS_EIGEN_MULTI_VECTOR_DOT_VECTOR_HPP_
#define CONTAINERS_SRC_OPS_EIGEN_MULTI_VECTOR_DOT_VECTOR_HPP_

namespace pressio{ namespace containers{ namespace ops{

/*
 * multi_vector dot vector
 */

// void specializiation for:
// * vec_type is a DYNAMIC eigen vector wrapper
// * result_vec_type is the same as vec_type
template <
  typename mvec_type,
  typename vec_type,
  typename result_vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_eigen<mvec_type>::value and
    containers::meta::is_vector_wrapper_eigen<vec_type>::value and
    containers::meta::is_vector_wrapper_eigen<result_vec_type>::value and
    containers::details::traits<result_vec_type>::is_dynamic
    > * = nullptr
  >
void dot(const mvec_type & mvA,
	 const vec_type & vecB,
	 result_vec_type & result)
{
  static_assert(containers::meta::are_scalar_compatible<mvec_type, vec_type, result_vec_type>::value,
		"Types are not scalar compatible");

  const auto numVecs = mvA.numVectors();
  // I can resize if needed because I know here it is a dynamic vector
  if ( result.extent(0) != numVecs )
    result.data()->resize(numVecs);
  *result.data() = (*mvA.data()).transpose() * (*vecB.data());
}


// non-void specialize for:
// * vec_type is a DYNAMIC eigen vector wrapper
// * result type is the same as vec_type
template <
  typename mvec_type,
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_eigen<mvec_type>::value and
    containers::meta::is_vector_wrapper_eigen<vec_type>::value and
    containers::details::traits<vec_type>::is_dynamic
    > * = nullptr
  >
vec_type dot(const mvec_type & mvA, const vec_type & vecB)
{
  static_assert(containers::meta::are_scalar_compatible<mvec_type, vec_type>::value,
		"Types are not scalar compatible");

  vec_type c(mvA.data()->cols());
  dot(mvA,vecB,c);
  return c;
}


// void specialize for:
// * vec_type is a generic eigen vector wrapper
// * result_vec_type is a STATIC eigen vector wrapper
template <
  typename mvec_type,
  typename vec_type,
  typename result_vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_eigen<mvec_type>::value and
    containers::meta::is_vector_wrapper_eigen<vec_type>::value and
    containers::meta::is_vector_wrapper_eigen<result_vec_type>::value and
    containers::details::traits<result_vec_type>::is_static 
    > * = nullptr
  >
void dot(const mvec_type & mvA, const vec_type & vecB, result_vec_type & result)
{

  static_assert(containers::meta::are_scalar_compatible<mvec_type, vec_type, result_vec_type>::value,
		"Types are not scalar compatible");

  // we are dealing with static vector type, so this needs to be true
  assert(result.extent(0) == mvA.data()->cols());
  // compute
  *result.data() = (*mvA.data()).transpose() * (*vecB.data());
}



// compute c += mvA^T vecB
// * vec_type is a DYNAMIC eigen vector wrapper
// * result_vec_type is an expression
template <
  typename mvec_type,
  typename vec_type,
  typename expr_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_eigen<mvec_type>::value and
    containers::meta::is_dynamic_vector_wrapper_eigen<vec_type>::value and
    containers::meta::are_scalar_compatible<mvec_type, vec_type, expr_type>::value and
    containers::meta::is_expression<expr_type>::value and
    ::pressio::containers::meta::is_vector_wrapper_eigen<
      typename ::pressio::containers::details::traits<expr_type>::data_t
      >::value
    > * = nullptr
  >
void updateWithDot(const mvec_type & mvA, const vec_type & vecB, expr_type & result)
{
  static_assert(containers::meta::are_scalar_compatible<mvec_type, vec_type, expr_type>::value,
		"Types are not scalar compatible");

  const auto numVecs = mvA.numVectors();
  assert(result.extent(0) == numVecs );
  result() += (*mvA.data()).transpose() * (*vecB.data());
}


}}}//end namespace pressio::ops
#endif
