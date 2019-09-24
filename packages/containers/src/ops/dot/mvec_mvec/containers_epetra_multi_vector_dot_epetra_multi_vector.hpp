/*
//@HEADER
// ************************************************************************
//
// containers_epetra_multi_vector_dot_epetra_multi_vector.hpp
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
#ifndef CONTAINERS_EPETRA_MULTI_VECTOR_DOT_EPETRA_MULTI_VECTOR_HPP_
#define CONTAINERS_EPETRA_MULTI_VECTOR_DOT_EPETRA_MULTI_VECTOR_HPP_

#include "../../containers_ops_meta.hpp"
#include "../../../multi_vector/containers_multi_vector_meta.hpp"

namespace pressio{ namespace containers{ namespace ops{

// Epetra multivector dot epetra multi vector

template <
  typename mvec_t,
  typename result_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_epetra<mvec_t>::value and
    ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<result_t>::value and
    ::pressio::containers::meta::wrapper_pair_have_same_scalar<mvec_t, result_t>::value
    > * = nullptr
  >
void dot(const mvec_t & mvA, const mvec_t & mvB, result_t & C){
  //using sc_t = typename containers::details::traits<mvec_t>::scalar_t;
  //  using eig_mat = Eigen::Matrix< sc_t, Eigen::Dynamic, Eigen::Dynamic>;
  //using res_t = containers::Matrix<eig_mat>;

  // how many vectors are in mvA and mvB
  const auto numVecsA = mvA.globalNumVectors();
  const auto numVecsB = mvB.globalNumVectors();
  assert( mvA.globalLength() == mvB.globalLength());
  auto const & mvAdata = *mvA.data();
  auto const & mvBdata = *mvB.data();

  assert(C.rows() == numVecsA);
  assert(C.cols() == numVecsB);
  // compute dot between every column of A with every col of B
  for (size_t i=0; i<(size_t)numVecsA; i++){
    for (size_t j=0; j<(size_t)numVecsB; j++){
      mvAdata(i)->Dot( *(mvBdata(j)), &C(i,j) );
    }
  }
}


template <
  typename mvec_t,
  typename result_t,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_epetra<mvec_t>::value and
    ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<result_t>::value and
    ::pressio::containers::meta::wrapper_pair_have_same_scalar<mvec_t, result_t>::value
    > * = nullptr
  >
result_t dot(const mvec_t & mvA, const mvec_t & mvB){
  // using sc_t = typename containers::details::traits<mvec_t>::scalar_t;
  // using eig_mat = Eigen::Matrix< sc_t, Eigen::Dynamic, Eigen::Dynamic>;
  // using res_t = containers::Matrix<eig_mat>;
  const auto numVecsA = mvA.globalNumVectors();
  const auto numVecsB = mvB.globalNumVectors();
  result_t C(numVecsA, numVecsB);
  dot(mvA, mvB, C);
  return C;
}


}}} // end namespace pressio::containers::ops
#endif
#endif
