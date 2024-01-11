/*
//@HEADER
// ************************************************************************
//
// ops_level2.hpp
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

#ifndef OPS_TPETRA_BLOCK_OPS_LEVEL2_HPP_
#define OPS_TPETRA_BLOCK_OPS_LEVEL2_HPP_

#include "Tpetra_idot.hpp"
#include <KokkosBlas1_axpby.hpp>
#include "KokkosBlas2_gemv.hpp"

namespace pressio{ namespace ops{

// -------------------------------
// y = beta * y + alpha*A*x
//
// x is Teuchos Vector
// A = tpetra block MultiVector
// y = tpetra block vector or column expression on mv
// -------------------------------
template < class A_type, class x_type, class y_type, class alpha_t, class beta_t>
std::enable_if_t<
  // level2 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<x_type>::rank == 1
  && ::pressio::Traits<y_type>::rank == 1
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra_block<A_type>::value
  && (   ::pressio::is_vector_tpetra_block<y_type>::value
      || ::pressio::is_expression_column_acting_on_tpetra_block<y_type>::value)
  && ::pressio::is_dense_vector_teuchos<x_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, x_type, y_type>::value
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  >
product(::pressio::nontranspose mode,
	const alpha_t & alpha,
	const A_type & A,
	const x_type & x,
	const beta_t & beta,
	y_type & y)
{

  auto A_tpetra_mv = A.getMultiVectorView();
  auto y_tpetra_v  = impl::get_underlying_tpetra_object(y);
  product(mode, alpha, A_tpetra_mv, x, beta, y_tpetra_v);
}

// -------------------------------
// y = beta * y + alpha*A*x
//
// x is Kokkos Vector
// A = tpetra block MultiVector
// y = tpetra block vector or column expression
// -------------------------------
template < class A_type, class x_type, class y_type, class alpha_t, class beta_t>
std::enable_if_t<
  // level2 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<x_type>::rank == 1
  && ::pressio::Traits<y_type>::rank == 1
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra_block<A_type>::value
  && (   ::pressio::is_vector_tpetra_block<y_type>::value
      || ::pressio::is_expression_column_acting_on_tpetra_block<y_type>::value)
  && (::pressio::is_vector_kokkos<x_type>::value
   || ::pressio::is_expression_acting_on_kokkos<x_type>::value)
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, x_type, y_type>::value
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  >
product(::pressio::nontranspose mode,
	const alpha_t & alpha,
	const A_type & A,
	const x_type & x,
	const beta_t & beta,
	y_type & y)
{

  auto A_tpetra_mv = A.getMultiVectorView();
  auto y_tpetra_v  = impl::get_underlying_tpetra_object(y);
  product(mode, alpha, A_tpetra_mv, x, beta, y_tpetra_v);
}

// -------------------------------
// y = beta * y + alpha*A*x
//
// x is Eigen Vector or Pressio expression based on Eigen container
// A = tpetra block MultiVector
// y = tpetra block vector or column expression
// -------------------------------
#ifdef PRESSIO_ENABLE_TPL_EIGEN
template < class A_type, class x_type, class y_type, class alpha_t, class beta_t>
std::enable_if_t<
  // level2 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<x_type>::rank == 1
  && ::pressio::Traits<y_type>::rank == 1
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra_block<A_type>::value
  && (   ::pressio::is_vector_tpetra_block<y_type>::value
      || ::pressio::is_expression_column_acting_on_tpetra_block<y_type>::value)
  && (::pressio::is_vector_eigen<x_type>::value
   || ::pressio::is_expression_acting_on_eigen<x_type>::value)
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, x_type, y_type>::value
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  >
product(::pressio::nontranspose mode,
	const alpha_t & alpha,
	const A_type & A,
	const x_type & x,
	const beta_t & beta,
	y_type & y)
{

  auto A_tpetra_mv = A.getMultiVectorView();
  auto y_tpetra_v  = impl::get_underlying_tpetra_object(y);
  product(mode, alpha, A_tpetra_mv, x, beta, y_tpetra_v);
}

// -------------------------------
// y = beta * y + alpha*A^T*x
//
// x = tpetra block Vector or column expression
// A = tpetra block MultiVector
// y = Eigen Vector or Pressio expression based on Eigen container
// -------------------------------
template < class A_type, class x_type, class y_type, class alpha_t, class beta_t>
std::enable_if_t<
  // level2 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<x_type>::rank == 1
  && ::pressio::Traits<y_type>::rank == 1
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra_block<A_type>::value
  && (   ::pressio::is_vector_tpetra_block<x_type>::value
      || ::pressio::is_expression_column_acting_on_tpetra_block<x_type>::value)
  && (::pressio::is_vector_eigen<y_type>::value
   || ::pressio::is_expression_acting_on_eigen<y_type>::value)
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, x_type, y_type>::value
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  >
product(::pressio::transpose mode,
	const alpha_t & alpha,
	const A_type & A,
	const x_type & x,
	const beta_t & beta,
	y_type & y)
{

  auto A_tpetra_mv = A.getMultiVectorView();
  auto x_tpetra_v  = impl::get_underlying_tpetra_object(x);
  product(mode, alpha, A_tpetra_mv, x_tpetra_v, beta, y);
}
#endif

// -------------------------------
// y = beta * y + alpha*A^T*x
//
// x = tpetra block Vector or column expression
// A = tpetra block MultiVector
// y = Kokkos vector
// -------------------------------
template < class A_type, class x_type, class y_type, class alpha_t, class beta_t>
std::enable_if_t<
  // level2 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<x_type>::rank == 1
  && ::pressio::Traits<y_type>::rank == 1
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra_block<A_type>::value
  && ( ::pressio::is_vector_tpetra_block<x_type>::value
    || ::pressio::is_expression_column_acting_on_tpetra_block<x_type>::value)
  && (::pressio::is_vector_kokkos<y_type>::value
   || ::pressio::is_expression_acting_on_kokkos<y_type>::value)
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, x_type, y_type>::value
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  >
product(::pressio::transpose mode,
	const alpha_t & alpha,
	const A_type & A,
	const x_type & x,
	const beta_t & beta,
	y_type & y)
{

  auto A_tpetra_mv = A.getMultiVectorView();
  auto x_tpetra_v  = impl::get_underlying_tpetra_object(x);
  product(mode, alpha, A_tpetra_mv, x_tpetra_v, beta, y);
}

}}//end namespace pressio::ops

#endif  // OPS_TPETRA_BLOCK_OPS_LEVEL2_HPP_
