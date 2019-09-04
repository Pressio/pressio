/*
//@HEADER
// ************************************************************************
//
// containers_eigen_ops_helper_impl.hpp
//                     		      Pressio 
// Copyright 2019 National Technology & Engineering Solutions of Sandia,LLC 
//							      (NTESS)
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

#ifndef CONTAINERS_EIGEN_OPS_HELPER_IMPL_HPP_
#define CONTAINERS_EIGEN_OPS_HELPER_IMPL_HPP_

#include "containers_ops_meta.hpp"
#include "../matrix/containers_matrix_meta.hpp"

namespace pressio{ namespace containers{ namespace ops{ namespace impl{


template <typename TA, typename TB, typename enable = void>
struct eigenMatMatProdRetTypeHelper;

// TA = dense, TB = dense, their product is dense
template <typename TA>
struct eigenMatMatProdRetTypeHelper<
  TA, TA,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_dense_matrix_wrapper_eigen<TA>::value
    >
  >{
  using prod_type = TA;
};

// TA = sparse, TB = sparse, their product is sparse
template <typename TA>
struct eigenMatMatProdRetTypeHelper<
  TA, TA,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_sparse_matrix_wrapper_eigen<TA>::value
    >
  >{
  using prod_type = TA;
};

// TA = sparse, TB = dense, their product is dense
template <typename TA, typename TB>
struct eigenMatMatProdRetTypeHelper<
  TA, TB,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_sparse_matrix_wrapper_eigen<TA>::value and
    containers::meta::is_dense_matrix_wrapper_eigen<TB>::value
    >
  >{
  using prod_type = TB;
};

// TA = dense, TB = sparse, their product is dense
template <typename TA, typename TB>
struct eigenMatMatProdRetTypeHelper<
  TA, TB,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_dense_matrix_wrapper_eigen<TA>::value and
    containers::meta::is_sparse_matrix_wrapper_eigen<TB>::value
    >
  >{
  using prod_type = TA;
};
//-------------------------------------------------------------


template <bool transpose_A, bool transpose_B>
struct eig_mat_mat_product;

// Default is C = A * B
template <>
struct eig_mat_mat_product<false, false>{

  template <typename TA, typename TB, typename TC>
  void operator()(const TA & A, const TB & B, TC & C) const{
    assert(A.data()->cols() == B.data()->rows());
    assert(C.data()->rows() == A.data()->rows());
    assert(C.data()->cols() == B.data()->cols());
    (*C.data()) = (*A.data()) * (*B.data());
  }

  template <typename TA, typename TB>
  auto operator()(const TA & A, const TB & B) const
    -> typename eigenMatMatProdRetTypeHelper<TA, TB>::prod_type
  {
    using prod_ret_t = typename eigenMatMatProdRetTypeHelper<TA, TB>::prod_type;
    prod_ret_t C(A.data()->rows(),B.data()->cols());
    (*this)(A,B,C);
    return C;
  }

};


// C = A^T * B
template <>
struct eig_mat_mat_product<true, false>{

  template <typename TA, typename TB, typename TC >
  void operator()(const TA & A, const TB & B, TC & C) const{
    assert(A.data()->rows() == B.data()->rows());
    assert(C.data()->rows() == A.data()->cols());
    assert(C.data()->cols() == B.data()->cols());
    (*C.data()) = (*A.data()).transpose() * (*B.data());
  }

  template <typename TA, typename TB>
  auto operator()(const TA & A, const TB & B) const
    -> typename eigenMatMatProdRetTypeHelper<TA, TB>::prod_type
  {
    using prod_ret_t = typename eigenMatMatProdRetTypeHelper<TA, TB>::prod_type;
    prod_ret_t C(A.data()->cols(),B.data()->cols());
    (*this)(A,B,C);
    return C;
  }
};


// C = A * B^T
template <>
struct eig_mat_mat_product<false, true>{

  template <typename TA, typename TB, typename TC >
  void operator()(const TA & A, const TB & B, TC & C) const{
    assert(A.data()->cols() == B.data()->cols());
    assert(C.data()->rows() == A.data()->rows());
    assert(C.data()->cols() == B.data()->rows());
    (*C.data()) = (*A.data()) * (*B.data()).transpose();
  }

  template <typename TA, typename TB>
  auto operator()(const TA & A, const TB & B) const
    -> typename eigenMatMatProdRetTypeHelper<TA, TB>::prod_type
  {
    using prod_ret_t = typename eigenMatMatProdRetTypeHelper<TA, TB>::prod_type;
    prod_ret_t C(A.data()->rows(), B.data()->rows());
    (*this)(A,B,C);
    return C;
  }
};


}}}}//end namespace pressio::containers::ops::impl
#endif
