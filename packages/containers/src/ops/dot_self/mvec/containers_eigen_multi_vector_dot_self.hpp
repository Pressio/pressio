/*
//@HEADER
// ************************************************************************
//
// containers_eigen_multi_vector_dot_self.hpp
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

#ifndef CONTAINERS_EIGEN_MULTI_VECTOR_DOT_SELF_HPP_
#define CONTAINERS_EIGEN_MULTI_VECTOR_DOT_SELF_HPP_

#include "../../containers_ops_meta.hpp"
#include "../../../multi_vector/containers_multi_vector_meta.hpp"

namespace pressio{ namespace containers{ namespace ops{

// Eigen multivector (A) dot self
// this is equivalent to doing A^T * A

template <
  typename mvec_t,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_eigen<mvec_t>::value
    > * = nullptr
  >
void dot_self(const mvec_t & A,
	      containers::Matrix<
	      Eigen::Matrix<typename containers::details::traits<mvec_t>::scalar_t,
	      Eigen::Dynamic, Eigen::Dynamic>
	      > & C)
{
  // // how many vectors are in A
  // auto numVecsA = A.numVectors();
  // auto const & Adata = *A.data();
  // assert(C.rows() == numVecsA);
  // assert(C.cols() == numVecsA);

  // // A dot A = A^T*A, which yields a symmetric matrix
  // // only need to compute half and fill remaining entries accordingly
  // for (auto i=0; i<numVecsA; i++){
  //   for (auto j=i; j<numVecsA; j++)
  //   {
  //     C(i,j) = A.data()->col(i).dot(A.data()->col(j));
  //     // fill the lower triangular part
  //     C(j,i) = C(i,j);
  //   }
  // }

  *C.data() = A.data()->transpose() * (*A.data());
}


template <
  typename mvec_t,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_eigen<mvec_t>::value
    > * = nullptr
  >
auto dot_self(const mvec_t & mvA)
  -> containers::Matrix<
    Eigen::Matrix<typename containers::details::traits<mvec_t>::scalar_t,
		  Eigen::Dynamic, Eigen::Dynamic>
    >{

  using sc_t = typename containers::details::traits<mvec_t>::scalar_t;
  using eig_mat = Eigen::Matrix< sc_t, Eigen::Dynamic, Eigen::Dynamic>;
  using res_t = containers::Matrix<eig_mat>;

  auto numVecsA = mvA.numVectors();
  res_t C(numVecsA, numVecsA);
  dot_self(mvA, C);
  return C;
}

}}} // end namespace pressio::containers::ops
#endif
