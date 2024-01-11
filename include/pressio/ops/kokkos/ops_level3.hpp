/*
//@HEADER
// ************************************************************************
//
// ops_level3.hpp
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

#ifndef OPS_KOKKOS_OPS_LEVEL3_HPP_
#define OPS_KOKKOS_OPS_LEVEL3_HPP_

#include "KokkosBlas3_gemm.hpp"

namespace pressio{ namespace ops{

/*
 * C = beta * C + alpha*op(A)*op(B)
*/

//-------------------------------------------
// specialize for op(A) = A and op(B) = B
//-------------------------------------------
template <typename A_type, typename B_type, typename alpha_type, typename beta_type, typename C_type>
std::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<B_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_dense_matrix_kokkos<A_type>::value
  && ::pressio::is_dense_matrix_kokkos<B_type>::value
  && ::pressio::is_dense_matrix_kokkos<C_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, B_type, C_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  && std::is_convertible<alpha_type, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_type,  typename ::pressio::Traits<A_type>::scalar_type>::value
  >
product(::pressio::nontranspose /*unused*/,
  ::pressio::nontranspose /*unused*/,
  const alpha_type alpha,
  const A_type & A,
  const B_type & B,
  const beta_type beta,
  C_type & C)
{

  assert( ::pressio::ops::extent(C, 0) == ::pressio::ops::extent(A, 0) );
  assert( ::pressio::ops::extent(C, 1) == ::pressio::ops::extent(B, 1) );
  assert( ::pressio::ops::extent(A, 1) == ::pressio::ops::extent(B, 0) );

  const char ctA = 'N';
  const char ctB = 'N';
  using sc_t = typename ::pressio::Traits<A_type>::scalar_type;
  const sc_t alpha_(alpha);
  const sc_t beta_(beta);
  ::KokkosBlas::gemm(&ctA, &ctB, alpha_, A, B, beta_, C);
}

//-------------------------------------------
// specialize for op(A) = A^T and op(B) = B
//-------------------------------------------
template <typename A_type, typename B_type, typename alpha_type, typename beta_type, typename C_type>
std::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<B_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_dense_matrix_kokkos<A_type>::value
  && ::pressio::is_dense_matrix_kokkos<B_type>::value
  && ::pressio::is_dense_matrix_kokkos<C_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, B_type, C_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  && std::is_convertible<alpha_type, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_type,  typename ::pressio::Traits<A_type>::scalar_type>::value
  >
product(::pressio::transpose /*unused*/,
	::pressio::nontranspose /*unused*/,
	const alpha_type alpha,
	const A_type & A,
	const B_type & B,
	const beta_type beta,
	C_type & C)
{

  assert( ::pressio::ops::extent(C, 0) == ::pressio::ops::extent(A, 1) );
  assert( ::pressio::ops::extent(C, 1) == ::pressio::ops::extent(B, 1) );
  assert( ::pressio::ops::extent(A, 0) == ::pressio::ops::extent(B, 0) );

  const char ctA = 'T';
  const char ctB = 'N';
  using sc_t = typename ::pressio::Traits<A_type>::scalar_type;
  const sc_t alpha_(alpha);
  const sc_t beta_(beta);
  ::KokkosBlas::gemm(&ctA, &ctB, alpha_, A, B, beta_, C);
}

/*----------------------------------------
* special case A==B and op(A)=A^T, op(B)=B
------------------------------------------*/
template <typename A_type, typename alpha_type, typename beta_type, typename C_type>
std::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_dense_matrix_kokkos<A_type>::value
  && ::pressio::is_dense_matrix_kokkos<C_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, C_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  && std::is_convertible<alpha_type, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_type,  typename ::pressio::Traits<A_type>::scalar_type>::value
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const alpha_type alpha,
	const A_type & A,
	const beta_type beta,
	C_type & C)
{
  product(modeA, modeB, alpha, A, A, beta, C);
}

/*---------------------------------------------------
* special case A==B and op(A)=A^T, op(B)=B, return C
-----------------------------------------------------*/
template <typename C_type, typename A_type, typename alpha_type>
std::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_dense_matrix_kokkos<A_type>::value
  && ::pressio::is_dense_matrix_kokkos<C_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, C_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  && std::is_convertible<alpha_type, typename ::pressio::Traits<A_type>::scalar_type>::value,
  C_type
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const alpha_type alpha,
	const A_type & A)
{
  using sc_t = typename ::pressio::Traits<A_type>::scalar_type;
  constexpr sc_t zero{0};
  C_type C("opsLev3C", A.extent(1), A.extent(1));
  product(modeA, modeB, alpha, A, A, zero, C);
  return C;
}

}}//end namespace pressio::ops
#endif  // OPS_KOKKOS_OPS_LEVEL3_HPP_
