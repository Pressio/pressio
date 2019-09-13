/*
//@HEADER
// ************************************************************************
//
// containers_eigen_multi_vector_product_eigen_vector.hpp
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

#ifndef CONTAINERS_CONTAINER_OPS_MVEC_VEC_PROD_EIGEN_MULTI_VECTOR_PRODUCT_EIGEN_VECTOR_HPP_
#define CONTAINERS_CONTAINER_OPS_MVEC_VEC_PROD_EIGEN_MULTI_VECTOR_PRODUCT_EIGEN_VECTOR_HPP_

#include "../containers_ops_meta.hpp"
#include "../../vector/containers_vector_meta.hpp"
#include "../../multi_vector/containers_multi_vector_meta.hpp"
#include "../../vector/concrete/containers_vector_sharedmem_eigen_dynamic.hpp"
#ifdef HAVE_TRILINOS
#include "../../vector/concrete/containers_vector_distributed_epetra.hpp"
#endif

namespace pressio{ namespace containers{ namespace ops{

//-----------------------------------------------------
// Eigen multivector product with eigen vector
template <typename mvec_type,
	  typename vec_type,
  ::pressio::mpl::enable_if_t<
   containers::meta::is_multi_vector_wrapper_eigen<mvec_type>::value and
   containers::meta::is_vector_wrapper_eigen<vec_type>::value and
   containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value
    > * = nullptr
  >
void product(const mvec_type & mvA,
	     const vec_type & vecB,
	     containers::Vector<Eigen::VectorXd> & C){

  assert( C.size() == mvA.length() );
  //zero out result
  C.setZero();
  // sizes of mvA
  auto numVecs = mvA.numVectors();
  auto Alength = mvA.length();

  // size of vecB
  assert(size_t(numVecs) == size_t(vecB.size()));

  // compute
  for (decltype(Alength) i=0; i<Alength; i++){
    for (decltype(numVecs) j=0; j<numVecs; j++){
      C[i] += mvA(i,j) * vecB[j];
    }
  }
}//end function


// result is constructed and returned
template <typename mvec_type,
	  typename vec_type,
  ::pressio::mpl::enable_if_t<
   containers::meta::is_multi_vector_wrapper_eigen<mvec_type>::value and
   containers::meta::is_vector_wrapper_eigen<vec_type>::value and
   containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value
    > * = nullptr
  >
auto product(const mvec_type & mvA, const vec_type & vecB)
-> containers::Vector<
    Eigen::Matrix<typename containers::details::traits<mvec_type>::scalar_t,
                  Eigen::Dynamic,1>
>{

  using sc_t = typename containers::details::traits<mvec_type>::scalar_t;
  containers::Vector<Eigen::Matrix<sc_t,Eigen::Dynamic,1>> c(mvA.length());
  product(mvA, vecB, c);
  return c;
}
//-------------------------------------------------------
//-------------------------------------------------------


}}}//end namespace pressio::containers::ops
#endif
