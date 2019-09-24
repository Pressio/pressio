/*
//@HEADER
// ************************************************************************
//
// containers_epetra_multi_vector_dot_epetra_vector.hpp
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

#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_EPETRA_MULTI_VECTOR_DOT_EPETRA_VECTOR_HPP_
#define CONTAINERS_EPETRA_MULTI_VECTOR_DOT_EPETRA_VECTOR_HPP_

#include "../../containers_ops_meta.hpp"
#include "../../../vector/containers_vector_meta.hpp"
#include "../../../multi_vector/containers_multi_vector_meta.hpp"

namespace pressio{ namespace containers{ namespace ops{

/*
 * overloads for c = A dot b
 * where:
 * A = wrapper of Epetra Multivector
 * b = epetra vector
 * c = a shared-mem vector, like eigen or armadillo
 */


//------------------------------------
// c = scalar *, passed in
//------------------------------------
template <
  typename mvec_type,
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_epetra<mvec_type>::value &&
    containers::meta::is_vector_wrapper_epetra<vec_type>::value &&
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value
    > * = nullptr
  >
void dot(const mvec_type & mvA,
	 const vec_type & vecB,
	 typename details::traits<mvec_type>::scalar_t * result){

  ///computes dot product of each vector in mvA
  ///with vecB storing each value in result
  /* Apparently, trilinos does not support this...
   * the native Dot() method of multivectors is only for
   * dot product of two multivectors of the same size.
   * So we have to extract each column vector
   * from mvA and do dot product one at a time
  */

  // how many vectors are in mvA
  const auto numVecs = mvA.globalNumVectors();
  auto * mvNatData = mvA.data();
  const auto * vecNatData = vecB.data();
  for (size_t i=0; i<(size_t)numVecs; i++){
    (*mvNatData)(i)->Dot(*vecNatData, &result[i]);
  }
}


//--------------------------------------------
// c = teuchos serial dense vector, passed in
//--------------------------------------------
template <
  typename mvec_type,
  typename vec_type,
  typename result_vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_epetra<mvec_type>::value and
    containers::meta::is_vector_wrapper_epetra<vec_type>::value and
    containers::meta::is_dense_vector_wrapper_teuchos<result_vec_type>::value and
    containers::meta::wrapper_triplet_have_same_scalar<mvec_type,
						       vec_type,
						       result_vec_type>::value and
    containers::details::traits<result_vec_type>::is_dynamic
    > * = nullptr
  >
void dot(const mvec_type & mvA,
	 const vec_type & vecB,
	 result_vec_type & result)
{
  const auto numVecs = mvA.globalNumVectors();
  if ( result.size() != numVecs )
    result.resize(numVecs);
  dot(mvA, vecB, result.data()->values());
}


//--------------------------------------
// c = Eigen DYNAMIC vector, passed in
//--------------------------------------
template <
  typename mvec_type,
  typename vec_type,
  typename result_vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_epetra<mvec_type>::value and
    containers::meta::is_vector_wrapper_epetra<vec_type>::value and
    containers::meta::is_vector_wrapper_eigen<result_vec_type>::value and
    containers::meta::wrapper_triplet_have_same_scalar<mvec_type,
						       vec_type,
						       result_vec_type>::value and
    containers::details::traits<result_vec_type>::is_dynamic
    > * = nullptr
  >
void dot(const mvec_type & mvA,
	 const vec_type & vecB,
	 result_vec_type & result)
{
  const auto numVecs = mvA.globalNumVectors();
  if ( result.size() != numVecs )
    result.resize(numVecs);
  dot(mvA, vecB, result.data()->data());
}


//--------------------------------------
// c = Eigen STATIC vector, passed in
//--------------------------------------
template <
  typename mvec_type,
  typename vec_type,
  typename result_vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_epetra<mvec_type>::value and
    containers::meta::is_vector_wrapper_epetra<vec_type>::value and
    containers::meta::is_vector_wrapper_eigen<result_vec_type>::value and
    containers::meta::wrapper_triplet_have_same_scalar<mvec_type,
						       vec_type,
						       result_vec_type>::value and
    containers::details::traits<result_vec_type>::is_static
    > * = nullptr
  >
void dot(const mvec_type & mvA,
	 const vec_type & vecB,
	 result_vec_type & result)
{
  assert(result.size() == mvA.globalNumVectors());
  dot(mvA, vecB, result.data()->data());
}


//--------------------------------------
// c = std::vector, passed in
//--------------------------------------
template <
  typename mvec_type,
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_epetra<mvec_type>::value &&
    containers::meta::is_vector_wrapper_epetra<vec_type>::value &&
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value
    > * = nullptr
  >
void dot(const mvec_type & mvA,
	 const vec_type & vecB,
	 std::vector<typename
	 details::traits<mvec_type>::scalar_t> & result){
  const auto numVecs = mvA.globalNumVectors();
  if ( result.size() != (size_t)numVecs )
    result.resize(numVecs);
  dot(mvA, vecB, result.data());
}


//--------------------------------------
// returns a std::vector
//--------------------------------------
template <
  typename mvec_type,
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_epetra<mvec_type>::value &&
    containers::meta::is_vector_wrapper_epetra<vec_type>::value &&
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value
    > * = nullptr
  >
std::vector<typename details::traits<mvec_type>::scalar_t>
dot(const mvec_type & mvA, const vec_type & vecB){

  using sc_t = typename details::traits<mvec_type>::scalar_t;
  // how many vectors are in mvA
  const auto numVecs = mvA.globalNumVectors();
  using res_t = std::vector<sc_t>;
  res_t res(numVecs);
  dot(mvA, vecB, res);
  return res;
}
//--------------------------------------------------------

}}} // end namespace pressio::containers::ops
#endif
#endif
