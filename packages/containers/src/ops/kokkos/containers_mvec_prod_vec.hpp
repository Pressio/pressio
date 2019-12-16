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

#include "../containers_ops_meta.hpp"
#include "../../multi_vector/containers_multi_vector_meta.hpp"
#include "KokkosBlas2_gemv.hpp"

namespace pressio{ namespace containers{ namespace ops{

/* multi_vector prod vector */

/* -------------------------------------------------------------------
 * specialize for kokkos mv wrapper operating on a kokkos vector wrapper
 *-------------------------------------------------------------------*/
template <
  typename mvec_type,
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_kokkos<mvec_type>::value and
    containers::meta::is_vector_wrapper_kokkos<vec_type>::value and
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value
    > * = nullptr
  >
void product(const mvec_type & A,
	 const vec_type & b,
	 containers::Vector<
	   Kokkos::View<
	     typename containers::details::traits<mvec_type>::scalar_t*,
	     typename containers::details::traits<mvec_type>::layout,
	     typename containers::details::traits<mvec_type>::execution_space
	   >
	 > & c)
{
  static_assert(meta::kokkos_wrapper_pair_have_same_exe_space<mvec_type, vec_type>::value,
		"product: MV and vec types need to have same execution space" );

  assert( A.data()->extent(1) == b.data()->extent(0) );
  assert( c.data()->extent(0) == A.data()->extent(0) );

  using sc_t = typename containers::details::traits<mvec_type>::scalar_t;
  constexpr auto zero = ::pressio::utils::constants::zero<sc_t>();
  constexpr auto one = ::pressio::utils::constants::one<sc_t>();

  const char ctA = 'N';
  KokkosBlas::gemv(&ctA, one, *A.data(), *b.data(), zero, *c.data());
}

template <
  typename mvec_type,
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_kokkos<mvec_type>::value and
    containers::meta::is_vector_wrapper_kokkos<vec_type>::value and
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value
    > * = nullptr
  >
auto product(const mvec_type & mvA, const vec_type & vecB)
-> containers::Vector<
    Kokkos::View<
      typename containers::details::traits<mvec_type>::scalar_t*,
      typename containers::details::traits<mvec_type>::layout,
      typename containers::details::traits<mvec_type>::execution_space
      >
  >{

  static_assert(meta::kokkos_wrapper_pair_have_same_exe_space<mvec_type, vec_type>::value,
		"product: MV and vec types need to have same execution space" );

  using sc_t = typename containers::details::traits<mvec_type>::scalar_t;
  using layout    = typename containers::details::traits<mvec_type>::layout;
  using exe_space = typename containers::details::traits<mvec_type>::execution_space;

  using v_t = Kokkos::View<sc_t*, layout, exe_space>;
  using res_t = containers::Vector<v_t>;

  res_t c("product_res", mvA.data()->extent(0));
  assert( mvA.data()->extent(1) == vecB.data()->extent(0) );
  product(mvA, vecB, c);
  return c;
}



/* -------------------------------------------------------------------
 * specialize for kokkos mv wrapper operating on col vector expression
 *-------------------------------------------------------------------*/
template <
  typename mvec_type,
  typename expr_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_kokkos<mvec_type>::value and
    ::pressio::containers::meta::is_expression<expr_type>::value
    > * = nullptr
  >
void product(const mvec_type & A,
	     const expr_type & exprObj,
	     containers::Vector<
	     Kokkos::View<
	     typename containers::details::traits<mvec_type>::scalar_t*,
	     typename containers::details::traits<mvec_type>::layout,
	     typename containers::details::traits<mvec_type>::execution_space
	     >
	     > & c)
{
  // type of data wrapped by the expression
  using expr_data_t = typename ::pressio::containers::details::traits<expr_type>::data_t;
  static_assert(meta::kokkos_wrapper_pair_have_same_exe_space<mvec_type, expr_data_t>::value,
  		"product: MV and expr types need to have same execution space" );

  assert( c.size() == A.length() );
  const auto numVecs = A.numVectors();
  assert(numVecs == exprObj.size());

  using sc_t = typename containers::details::traits<mvec_type>::scalar_t;
  constexpr auto zero = ::pressio::utils::constants::zero<sc_t>();
  constexpr auto one = ::pressio::utils::constants::one<sc_t>();
  const char ctA = 'N';
  KokkosBlas::gemv(&ctA, one, *A.data(), exprObj(), zero, *c.data());
}

template <
  typename mvec_type,
  typename expr_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_kokkos<mvec_type>::value and
    ::pressio::containers::meta::is_expression<expr_type>::value
    > * = nullptr
  >
auto product(const mvec_type & mvA, const expr_type & exprObj)
-> containers::Vector<
    Kokkos::View<
      typename containers::details::traits<mvec_type>::scalar_t*,
      typename containers::details::traits<mvec_type>::layout,
      typename containers::details::traits<mvec_type>::execution_space
      >
  >
{
  // type of data wrapped by the expression
  using expr_data_t = typename ::pressio::containers::details::traits<expr_type>::data_t;
  static_assert(meta::kokkos_wrapper_pair_have_same_exe_space<mvec_type, expr_data_t>::value,
  		"product: MV and expr types need to have same execution space" );

  using sc_t = typename containers::details::traits<mvec_type>::scalar_t;
  using layout    = typename containers::details::traits<mvec_type>::layout;
  using exe_space = typename containers::details::traits<mvec_type>::execution_space;

  using v_t = Kokkos::View<sc_t*, layout, exe_space>;
  using res_t = containers::Vector<v_t>;

  res_t c("product_res", mvA.data()->extent(0));
  product(mvA, exprObj, c);
  return c;
}


}}}//end namespace pressio::containers::ops
#endif
#endif
