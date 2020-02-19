/*
//@HEADER
// ************************************************************************
//
// containers_mvec_prod_mvec.hpp
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

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#ifndef CONTAINERS_SRC_OPS_KOKKOS_MULTI_VECTOR_PROD_MVEC_HPP_
#define CONTAINERS_SRC_OPS_KOKKOS_MULTI_VECTOR_PROD_MVEC_HPP_

#include "KokkosBlas3_gemm.hpp"

namespace pressio{ namespace containers{ namespace ops{

/*
 * for tpetra:
 *
 * C = beta * C + alpha*op(A)*op(B)
 *
*/

template <typename A_type, typename B_type, typename scalar_type, typename C_type>
::pressio::mpl::enable_if_t<
  ::pressio::containers::meta::is_multi_vector_wrapper_kokkos<A_type>::value and
  ::pressio::containers::meta::is_multi_vector_wrapper_kokkos<B_type>::value and
  ::pressio::containers::meta::is_matrix_wrapper_kokkos<C_type>::value
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const A_type & A,
	const B_type & B,
	const scalar_type beta,
	C_type & C)
{
  static_assert(containers::meta::are_scalar_compatible<A_type, B_type, C_type>::value,
		"Types are not scalar compatible");

  // // std::is_same<
  // //   typename containers::details::traits<result_t>::execution_space,
  // //   typename containers::details::traits<mvec_t>::execution_space
  // // >::value and
  // // std::is_same<
  // //   typename containers::details::traits<result_t>::layout,
  // //   typename containers::details::traits<mvec_t>::layout
  // // >::value

  const char ctA = 'T';
  const char ctB = 'N';
  KokkosBlas::gemm(&ctA, &ctB, alpha, *A.data(), *B.data(), beta, *C.data());
}


/***********************************
* special case A==B (for now use impl above)
**********************************/
template <typename A_type, typename scalar_type, typename C_type>
::pressio::mpl::enable_if_t<
  ::pressio::containers::meta::is_multi_vector_wrapper_kokkos<A_type>::value and
  ::pressio::containers::meta::is_matrix_wrapper_kokkos<C_type>::value
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const A_type & A,
	const scalar_type beta,
	C_type & C)
{
  product(modeA, modeB, alpha, A, A, beta, C);
}

template <typename C_type, typename A_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  ::pressio::containers::meta::is_multi_vector_wrapper_kokkos<A_type>::value and
  ::pressio::containers::meta::is_matrix_wrapper_kokkos<C_type>::value,
  C_type
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const A_type & A)
{
  constexpr auto zero = ::pressio::utils::constants::zero<scalar_type>();
  C_type C(A.numVectors(), A.numVectors());
  product(modeA, modeB, alpha, A, A, zero, C);
  return C;
}

}}}//end namespace pressio::containers::ops
#endif
#endif






// // C += A^T B
// template <
//   typename mvec_t,
//   typename expr_t,
//   ::pressio::mpl::enable_if_t<
//     containers::meta::is_multi_vector_wrapper_kokkos<mvec_t>::value and
//     containers::meta::is_expression<expr_t>::value and
//     ::pressio::containers::meta::is_matrix_wrapper_kokkos<
//       typename ::pressio::containers::details::traits<expr_t>::data_t
//       >::value and
//     std::is_same<
//       typename containers::details::traits<expr_t>::execution_space,
//       typename containers::details::traits<mvec_t>::execution_space
//     >::value
//     > * = nullptr
//   >
// void updateWithDot(const mvec_t & A, const mvec_t & B, expr_t & C)
// {
//   static_assert(containers::meta::are_scalar_compatible<mvec_t, expr_t>::value,
//     "Types are not scalar compatible");

//   using sc_t = typename containers::details::traits<mvec_t>::scalar_t;
//   constexpr auto zero = ::pressio::utils::constants::zero<sc_t>();
//   constexpr auto one = ::pressio::utils::constants::one<sc_t>();
//   const char ctA = 'T';
//   const char ctB = 'N';
//   KokkosBlas::gemm(&ctA, &ctB, one, *A.data(), *B.data(), one, C());
// }

// /* ----------------------------------------------
//  * result_t = an expression
//  ---------------------------------------------- */
// template <
//   typename mvec_t,
//   typename result_t,
//   ::pressio::mpl::enable_if_t<
//     containers::meta::is_multi_vector_wrapper_kokkos<mvec_t>::value and
//     containers::meta::is_expression<result_t>::value
//     > * = nullptr
//   >
// void dot(const mvec_t & A, const mvec_t & B, result_t & C)
// {
//   static_assert(containers::meta::are_scalar_compatible<mvec_t, result_t>::value,
//     "Types are not scalar compatible");

//   using sc_t = typename containers::details::traits<mvec_t>::scalar_t;
//   constexpr auto zero = ::pressio::utils::constants::zero<sc_t>();
//   constexpr auto one = ::pressio::utils::constants::one<sc_t>();
//   const char ctA = 'T';
//   const char ctB = 'N';
//   KokkosBlas::gemm( &ctA, &ctB, one, *A.data(), *B.data(), zero, C() );
// }
