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

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#ifndef CONTAINERS_SRC_OPS_EPETRA_MULTI_VECTOR_DOT_VECTOR_HPP_
#define CONTAINERS_SRC_OPS_EPETRA_MULTI_VECTOR_DOT_VECTOR_HPP_

namespace pressio{ namespace containers{ namespace ops{

/*
 * multi_vector dot vector
 *
 * overloads for c = A dot b
 * where:
 * A = wrapper of Epetra Multivector
 * b = epetra vector
 * c = a shared-mem vector, like eigen or armadillo or std::vector
 */


//------------------------------------
// c = scalar *, passed in
//------------------------------------
template <
  typename mvec_type,
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_epetra<mvec_type>::value &&
    containers::meta::is_vector_wrapper_epetra<vec_type>::value
    > * = nullptr
  >
void dot(const mvec_type & mvA,
	 const vec_type & vecB,
	 typename details::traits<mvec_type>::scalar_t * result)
{
  static_assert(containers::meta::wrappers_have_same_scalar<mvec_type, vec_type>::value,
    "Types are not scalar compatible");

  ///computes dot product of each vector in mvA
  ///with vecB storing each value in result
  /* Apparently, trilinos does not support this...
   * the native Dot() method of multivectors is only for
   * dot product of two multivectors of the same size.
   * So we have to extract each column vector
   * from mvA and do dot product one at a time
  */

  // how many vectors are in mvA
  const auto numVecs = mvA.numVectorsGlobal();
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
    containers::details::traits<result_vec_type>::is_dynamic
    > * = nullptr
  >
void dot(const mvec_type & mvA,
	 const vec_type & vecB,
	 result_vec_type & result)
{
  static_assert(containers::meta::wrappers_have_same_scalar<mvec_type, vec_type, result_vec_type>::value,
    "Types are not scalar compatible");

  const auto numVecs = mvA.numVectorsGlobal();
  result.data()->resize(numVecs);
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
    containers::meta::is_dynamic_vector_wrapper_eigen<result_vec_type>::value
    > * = nullptr
  >
void dot(const mvec_type & mvA,
	 const vec_type & vecB,
	 result_vec_type & result)
{
  static_assert(containers::meta::wrappers_have_same_scalar<mvec_type, vec_type, result_vec_type>::value,
    "Types are not scalar compatible");

  const auto numVecs = mvA.numVectorsGlobal();
  if ( result.extent(0) != numVecs ){
    result.data()->resize(numVecs);
  }
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
    containers::details::traits<result_vec_type>::is_static
    > * = nullptr
  >
void dot(const mvec_type & mvA,
	 const vec_type & vecB,
	 result_vec_type & result)
{
  static_assert(containers::meta::wrappers_have_same_scalar<mvec_type, vec_type, result_vec_type>::value,
    "Types are not scalar compatible");

  assert(result.extent(0) == mvA.numVectorsGlobal());
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
    containers::meta::is_vector_wrapper_epetra<vec_type>::value
    > * = nullptr
  >
void dot(const mvec_type & mvA,
	 const vec_type & vecB,
	 std::vector<typename
	 details::traits<mvec_type>::scalar_t> & result)
{
  static_assert(containers::meta::wrappers_have_same_scalar<mvec_type, vec_type>::value,
    "Types are not scalar compatible");

  const auto numVecs = mvA.numVectorsGlobal();
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
    containers::meta::is_vector_wrapper_epetra<vec_type>::value
    > * = nullptr
  >
std::vector<typename details::traits<mvec_type>::scalar_t>
dot(const mvec_type & mvA, const vec_type & vecB)
{
  static_assert(containers::meta::wrappers_have_same_scalar<mvec_type, vec_type>::value,
    "Types are not scalar compatible");

  using sc_t = typename details::traits<mvec_type>::scalar_t;
  const auto numVecs = mvA.numVectorsGlobal();
  using res_t = std::vector<sc_t>;

  res_t res(numVecs);
  dot(mvA, vecB, res);
  return res;
}
//--------------------------------------------------------

}}}//end namespace pressio::containers::ops
#endif
#endif
