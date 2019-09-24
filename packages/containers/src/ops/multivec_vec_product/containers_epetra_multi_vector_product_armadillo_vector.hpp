/*
//@HEADER
// ************************************************************************
//
// containers_epetra_multi_vector_product_armadillo_vector.hpp
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

#if defined(HAVE_TRILINOS) && defined(HAVE_ARMADILLO)
#ifndef CONTAINERS_CONTAINER_OPS_MVEC_VEC_PROD_EPETRA_MULTI_VECTOR_PRODUCT_ARMADILLO_VECTOR_HPP_
#define CONTAINERS_CONTAINER_OPS_MVEC_VEC_PROD_EPETRA_MULTI_VECTOR_PRODUCT_ARMADILLO_VECTOR_HPP_

#include "../containers_ops_meta.hpp"
#include "../../vector/containers_vector_meta.hpp"
#include "../../multi_vector/containers_multi_vector_meta.hpp"
#include "../../vector/concrete/containers_vector_sharedmem_eigen_dynamic.hpp"
#include "../../vector/concrete/containers_vector_distributed_epetra.hpp"

namespace pressio{ namespace containers{ namespace ops{

//-----------------------------------------------------
//  Epetra multivector with eigen or armadillo vector
// we pass the result object
template <
  typename mvec_type,
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_epetra<mvec_type>::value and
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value and
    (containers::meta::is_column_vector_wrapper_armadillo<vec_type>::value or
     containers::meta::is_row_vector_wrapper_armadillo<vec_type>::value)
    > * = nullptr
  >
void product(const mvec_type & mvA,
	     const vec_type & vecB,
	     containers::Vector<Epetra_Vector> & C){

  //zero out result
  C.setZero();
  // how many vectors are in mvA
  const auto numVecs = mvA.globalNumVectors();
  // size of vecB
  const size_t vecBLen = vecB.size();
  assert(size_t(numVecs) == vecBLen);
  // the data map of the multivector
  const auto mvMap = mvA.getDataMap();
  // my number of rows
  const auto myNrows = mvMap.NumMyElements();

  // loop
  for (decltype(myNrows) i=0; i<myNrows; i++){
    for (decltype(numVecs) j=0; j<numVecs; j++){
      C[i] += mvA(i,j) * vecB[j];
    }
  }
}
//-------------------------------------------------------

// result is returned
template <typename mvec_type,
	  typename vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_epetra<mvec_type>::value and
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value and
    (containers::meta::is_column_vector_wrapper_armadillo<vec_type>::value or
     containers::meta::is_row_vector_wrapper_armadillo<vec_type>::value)
  > * = nullptr
 >
containers::Vector<Epetra_Vector>
product(const mvec_type & mvA, const vec_type & vecB) {

  // here, mvA is distrubted, but vecB is NOT.
  // we interpret this as a linear combination of vectors

  // the data map of the multivector
  const auto mvMap = mvA.getDataMap();
  // result is an Epetra Vector with same distribution of mvA
  using res_t = containers::Vector<Epetra_Vector>;
  res_t c(mvMap);
  product(mvA, vecB, c);
  return c;
}
//-------------------------------------------------------

}}}//end namespace pressio::containers::ops
#endif
#endif //HAVE_TRILINOS && HAVE_ARMADILLO
