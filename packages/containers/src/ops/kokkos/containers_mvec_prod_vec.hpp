/*
//@HEADER
// ************************************************************************
//
// containers_mvec_prod_vec.hpp
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
#ifndef CONTAINERS_SRC_OPS_KOKKOS_MULTI_VECTOR_PROD_VECTOR_HPP_
#define CONTAINERS_SRC_OPS_KOKKOS_MULTI_VECTOR_PROD_VECTOR_HPP_

#include "KokkosBlas2_gemv.hpp"

namespace pressio{ namespace containers{ namespace ops{

/*
 * multi_vector prod vector
 *
 * y = beta * y + alpha*op(A)*x
 *
*/

template < typename A_type, typename x_type, typename scalar_type, typename y_type>
::pressio::mpl::enable_if_t<
  containers::meta::is_multi_vector_wrapper_kokkos<A_type>::value and
  containers::meta::is_vector_wrapper_kokkos<x_type>::value
  >
product(::pressio::nontranspose mode,
	const scalar_type alpha,
	const A_type & A,
	const x_type & x,
	const scalar_type beta,
	y_type & y)
{
  static_assert(containers::meta::are_scalar_compatible<A_type, x_type, y_type>::value,
		"Types are not scalar compatible");
  static_assert(meta::have_matching_execution_space<A_type, x_type, y_type>::value,
		"operands need to have same execution space" );

  assert( y.data()->extent(0) == A.data()->extent(0) );
  assert( A.data()->extent(1) == x.data()->extent(0) );
  const char ctA = 'N';
  KokkosBlas::gemv( &ctA, alpha, *A.data(), *x.data(), beta, *y.data() );
}


template < typename A_type, typename x_type, typename scalar_type, typename y_type>
::pressio::mpl::enable_if_t<
  containers::meta::is_multi_vector_wrapper_kokkos<A_type>::value and
  containers::meta::is_vector_wrapper_kokkos<x_type>::value
  >
product(::pressio::transpose mode,
	const scalar_type alpha,
	const A_type & A,
	const x_type & x,
	const scalar_type beta,
	y_type & y)
{
  static_assert(containers::meta::are_scalar_compatible<A_type, x_type, y_type>::value,
		"Types are not scalar compatible");
  static_assert(meta::have_matching_execution_space<A_type, x_type, y_type>::value,
		"operands need to have same execution space" );

  assert( y.data()->extent(0) == A.data()->extent(1) );
  assert( A.data()->extent(0) == x.data()->extent(0) );
  const char ctA = 'T';
  KokkosBlas::gemv( &ctA, alpha, *A.data(), *x.data(), beta, *y.data() );
}

}}}//end namespace pressio::containers::ops
#endif
#endif


// /* -------------------------------------------------------------------
//  * specialize for kokkos mv wrapper operating on col vector expression
//  *-------------------------------------------------------------------*/
// template <
//   typename mvec_type,
//   typename expr_type,
//   ::pressio::mpl::enable_if_t<
//     containers::meta::is_multi_vector_wrapper_kokkos<mvec_type>::value and
//     ::pressio::containers::meta::is_expression<expr_type>::value
//     > * = nullptr
//   >
// void product(const mvec_type & A,
// 	     const expr_type & exprObj,
// 	     containers::Vector<
// 	     Kokkos::View<
// 	     typename containers::details::traits<mvec_type>::scalar_t*,
// 	     typename containers::details::traits<mvec_type>::layout,
// 	     typename containers::details::traits<mvec_type>::execution_space
// 	     >
// 	     > & c)
// {
//   static_assert(containers::meta::are_scalar_compatible<mvec_type, expr_type>::value,
//     "Types are not scalar compatible");

//   // type of data wrapped by the expression
//   using expr_data_t = typename ::pressio::containers::details::traits<expr_type>::data_t;
//   static_assert(meta::have_matching_execution_space<mvec_type, expr_data_t>::value,
//   		"product: MV and expr types need to have same execution space" );

//   assert( c.extent(0) == A.extent(0) );
//   const auto numVecs = A.numVectors();
//   assert(numVecs == exprObj.extent(0));

//   using sc_t = typename containers::details::traits<mvec_type>::scalar_t;
//   constexpr auto zero = ::pressio::utils::constants::zero<sc_t>();
//   constexpr auto one = ::pressio::utils::constants::one<sc_t>();
//   const char ctA = 'N';
//   KokkosBlas::gemv(&ctA, one, *A.data(), exprObj(), zero, *c.data());
// }

// template <
//   typename mvec_type,
//   typename expr_type,
//   ::pressio::mpl::enable_if_t<
//     containers::meta::is_multi_vector_wrapper_kokkos<mvec_type>::value and
//     ::pressio::containers::meta::is_expression<expr_type>::value
//     > * = nullptr
//   >
// auto product(const mvec_type & mvA, const expr_type & exprObj)
// -> containers::Vector<
//     Kokkos::View<
//       typename containers::details::traits<mvec_type>::scalar_t*,
//       typename containers::details::traits<mvec_type>::layout,
//       typename containers::details::traits<mvec_type>::execution_space
//       >
//   >
// {
//   static_assert(containers::meta::are_scalar_compatible<mvec_type, expr_type>::value,
//     "Types are not scalar compatible");

//   // type of data wrapped by the expression
//   using expr_data_t = typename ::pressio::containers::details::traits<expr_type>::data_t;
//   static_assert(meta::have_matching_execution_space<mvec_type, expr_data_t>::value,
//   		"product: MV and expr types need to have same execution space" );

//   using sc_t = typename containers::details::traits<mvec_type>::scalar_t;
//   using layout    = typename containers::details::traits<mvec_type>::layout;
//   using exe_space = typename containers::details::traits<mvec_type>::execution_space;

//   using v_t = Kokkos::View<sc_t*, layout, exe_space>;
//   using res_t = containers::Vector<v_t>;

//   res_t c("product_res", mvA.data()->extent(0));
//   product(mvA, exprObj, c);
//   return c;
// }
