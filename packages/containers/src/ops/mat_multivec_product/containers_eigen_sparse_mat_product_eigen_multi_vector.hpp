/*
//@HEADER
// ************************************************************************
//
// containers_eigen_sparse_mat_product_eigen_multi_vector.hpp
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

#ifndef CONTAINERS_CONTAINERS_OPS_EIGEN_SPARSE_MAT_PRODUCT_EIGEN_MULTI_VECTOR_HPP_
#define CONTAINERS_CONTAINERS_OPS_EIGEN_SPARSE_MAT_PRODUCT_EIGEN_MULTI_VECTOR_HPP_

#include "../containers_ops_meta.hpp"
#include "../../vector/containers_vector_meta.hpp"
#include "../../multi_vector/containers_multi_vector_meta.hpp"
#include "../../matrix/containers_matrix_meta.hpp"
#include "../../matrix/concrete/containers_matrix_sparse_sharedmem_eigen.hpp"
#include "../../matrix/concrete/containers_matrix_dense_sharedmem_eigen_dynamic.hpp"
#include "../../multi_vector/concrete/containers_multi_vector_sharedmem_eigen_dynamic.hpp"

namespace pressio{ namespace containers{ namespace ops{

/*---------------------------------------------
 * C = A * B
 * A : eigen sparse matrix wrapper
 * B: eigen multivector wrapper
 * C: eigen dense matrix wrapper
 *-----------------------------------------------*/
template <typename mat_type,
	  typename mvec_type,
  ::pressio::mpl::enable_if_t<
   ::pressio::containers::meta::is_sparse_matrix_wrapper_eigen<mat_type>::value and
   ::pressio::containers::meta::is_multi_vector_wrapper_eigen<mvec_type>::value and
    ::pressio::containers::meta::wrapper_pair_have_same_scalar<mat_type, mvec_type>::value
    > * = nullptr
  >
void product(const mat_type & A,
	     const mvec_type & mv,
	     ::pressio::containers::Matrix<Eigen::MatrixXd> & C){

  assert( C.rows() == A.rows() );
  assert( mv.length() == A.cols() );
  assert( C.cols() == mv.numVectors() );
  (*C.data()) = (*A.data()) * (*mv.data());
}//end function

// construct and return result
template <typename mat_type,
	  typename mvec_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_sparse_matrix_wrapper_eigen<mat_type>::value and
    ::pressio::containers::meta::is_multi_vector_wrapper_eigen<mvec_type>::value and
    ::pressio::containers::meta::wrapper_pair_have_same_scalar<mat_type, mvec_type>::value
    > * = nullptr
  >
auto product(const mat_type & A, const mvec_type & mv)
  -> ::pressio::containers::Matrix<Eigen::MatrixXd>
{
  ::pressio::containers::Matrix<Eigen::MatrixXd> C(A.rows(), mv.numVectors());
  product(A,mv,C);
  return C;
}//end function


}}}//end namespace pressio::containers::ops
#endif
