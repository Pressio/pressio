/*
//@HEADER
// ************************************************************************
//
// containers_eigen_mat_product_mat.hpp
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

#ifndef CONTAINERS_EIGEN_MAT_PRODUCT_MAT_HPP_
#define CONTAINERS_EIGEN_MAT_PRODUCT_MAT_HPP_

#include "../containers_ops_meta.hpp"
#include "../../matrix/containers_matrix_meta.hpp"
#include "../containers_eigen_ops_helper_impl.hpp"

namespace pressio{ namespace containers{ namespace ops{

/*---------------------------------------------------------
 *
 * C = A B
 * A matrix wrapper for an eigen matrix (static/dynamic)
 * B matrix wrapper for an eigen matrix (static/dynamic)
 *
 *--------------------------------------------------------*/


template <
  typename TA, typename TB, typename TC,
  bool transposeA = false, bool transposeB = false,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_matrix_wrapper_eigen<TA>::value and
    containers::meta::is_matrix_wrapper_eigen<TB>::value and
    containers::meta::is_matrix_wrapper_eigen<TC>::value and
    containers::meta::wrapper_triplet_have_same_scalar<TA,TB,TC>::value
    > * = nullptr
  >
void product(const TA & A, const TB & B, TC & C){

  using implClass_t = impl::eig_mat_mat_product<transposeA,
						transposeB>;
  implClass_t()(A,B,C);
}

// DENSE times DENSE
template <
  typename T1, typename T2,
    bool transposeA = false, bool transposeB = false,
    ::pressio::mpl::enable_if_t<
      containers::meta::is_dense_matrix_wrapper_eigen<T1>::value and
      containers::meta::is_dense_matrix_wrapper_eigen<T2>::value and
      containers::meta::wrapper_pair_have_same_scalar<T1,T2>::value
      > * = nullptr
    >
auto product(const T1 & A, const T2 & B)
  -> typename impl::eigenMatMatProdRetTypeHelper<T1, T2>::prod_type
{
  using implClass_t = impl::eig_mat_mat_product<transposeA,
						transposeB>;
  return implClass_t()(A,B);
}

// SPARSE times DENSE
template <
  typename T1, typename T2,
    bool transposeA = false, bool transposeB = false,
    ::pressio::mpl::enable_if_t<
      containers::meta::is_sparse_matrix_wrapper_eigen<T1>::value and
      containers::meta::is_dense_matrix_wrapper_eigen<T2>::value and
      containers::meta::wrapper_pair_have_same_scalar<T1,T2>::value
      > * = nullptr
    >
auto product(const T1 & A, const T2 & B)
  -> typename impl::eigenMatMatProdRetTypeHelper<T1, T2>::prod_type
{
  using implClass_t = impl::eig_mat_mat_product<transposeA,
						transposeB>;
  return implClass_t()(A,B);
}

// DENSE times SPARSE
template <
  typename T1, typename T2,
    bool transposeA = false, bool transposeB = false,
    ::pressio::mpl::enable_if_t<
      containers::meta::is_dense_matrix_wrapper_eigen<T1>::value and
      containers::meta::is_sparse_matrix_wrapper_eigen<T2>::value and
      containers::meta::wrapper_pair_have_same_scalar<T1,T2>::value
      > * = nullptr
    >
auto product(const T1 & A, const T2 & B)
  -> typename impl::eigenMatMatProdRetTypeHelper<T1, T2>::prod_type
{
  using implClass_t = impl::eig_mat_mat_product<transposeA,
						transposeB>;
  return implClass_t()(A,B);
}


// SPARSE times SPARSE
template <
  typename T1, typename T2,
    bool transposeA = false, bool transposeB = false,
    ::pressio::mpl::enable_if_t<
      containers::meta::is_sparse_matrix_wrapper_eigen<T1>::value and
      containers::meta::is_sparse_matrix_wrapper_eigen<T2>::value and
      containers::meta::wrapper_pair_have_same_scalar<T1,T2>::value
      > * = nullptr
    >
auto product(const T1 & A, const T2 & B)
  -> typename impl::eigenMatMatProdRetTypeHelper<T1, T2>::prod_type
{
  using implClass_t = impl::eig_mat_mat_product<transposeA,
						transposeB>;
  return implClass_t()(A,B);
}


}}} // end namespace pressio::containers::ops
#endif
