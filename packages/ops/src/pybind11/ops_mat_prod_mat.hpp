/*
//@HEADER
// ************************************************************************
//
// ops_mat_prod_mat.hpp
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

#ifndef OPS_PYBIND11_OPS_MAT_PROD_MAT_HPP_
#define OPS_PYBIND11_OPS_MAT_PROD_MAT_HPP_

namespace pressio{ namespace ops{

/*
 * C = beta * C + alpha*op(A)*op(B)
 *
*/

/***********************************
* special case A==B and op(A) = transpose
* i.e.: C = beta * C + alpha*A^T*A
**********************************/

template <typename A_type, typename scalar_type, typename C_type>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_fstyle_matrix_wrapper_pybind<A_type>::value and
  ::pressio::containers::predicates::is_fstyle_matrix_wrapper_pybind<C_type>::value
  >
product(::pressio::transpose modeA,
  ::pressio::nontranspose modeB,
  const scalar_type alpha,
  const A_type & A,
  const scalar_type beta,
  C_type & C)
{
  static_assert(containers::predicates::are_scalar_compatible<A_type, C_type>::value,
    "Types are not scalar compatible");

  // // we know that operators passed in are f-style arrays, so use eigenmap to manage them
  // // and perform computations (this should be cheaper than using scipy blas below)
  // using mat_t = Eigen::Matrix<scalar_type, -1, -1, Eigen::ColMajor>;
  // Eigen::Map<const mat_t> Am(A.data()->data(), A.extent(0), A.extent(1));
  // auto AmT = Am.transpose();
  // Eigen::Map<mat_t> Cm(C.data()->mutable_data(), C.extent(0), C.extent(1));
  // Cm = beta * Cm + alpha * AmT * Am;

  // NOTE: need to check if doing this import is expensive,
  // and assess whether we can use blas directly when we know
  // that objects involved are dense with not strange layout.
  pybind11::object pyblas = pybind11::module::import("scipy.linalg.blas");
  constexpr auto one  = ::pressio::utils::constants<scalar_type>::one();
  constexpr auto no   = ::pressio::utils::constants<int>::zero();
  constexpr auto yes  = ::pressio::utils::constants<int>::one();
  constexpr auto transA = yes;
  constexpr auto transB = no;
  constexpr auto ovw    = yes;
  pyblas.attr("dgemm")(one, *A.data(), *A.data(), beta, *C.data(), transA, transB, ovw);
}

template <typename C_type, typename A_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_fstyle_matrix_wrapper_pybind<A_type>::value and
  ::pressio::containers::predicates::is_fstyle_matrix_wrapper_pybind<C_type>::value,
  C_type
  >
product(::pressio::transpose modeA,
  ::pressio::nontranspose modeB,
  const scalar_type alpha,
  const A_type & A)
{
  static_assert(containers::predicates::are_scalar_compatible<A_type, C_type>::value,
    "Types are not scalar compatible");

  constexpr auto zero  = ::pressio::utils::constants<scalar_type>::zero();
  C_type C(A.extent(1), A.extent(1));
  product(modeA, modeB, alpha, A, zero, C);
  return C;
}


}}//end namespace pressio::ops
#endif  // OPS_PYBIND11_OPS_MAT_PROD_MAT_HPP_
