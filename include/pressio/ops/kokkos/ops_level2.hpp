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

#ifndef OPS_KOKKOS_OPS_LEVEL2_HPP_
#define OPS_KOKKOS_OPS_LEVEL2_HPP_

#include "KokkosBlas2_gemv.hpp"

namespace pressio{ namespace ops{

/*
 * multi_vector prod vector
 *
 * y = beta * y + alpha*op(A)*x
 *
*/

//-------------------------------
// specialize for op(A) = A
//-------------------------------
template < typename A_type, typename x_type, typename scalar_type, typename y_type>
::pressio::mpl::enable_if_t<
  (::pressio::is_native_container_kokkos<A_type>::value or
   ::pressio::is_expression_acting_on_kokkos<A_type>::value) and
  (::pressio::is_native_container_kokkos<x_type>::value or
   ::pressio::is_expression_acting_on_kokkos<x_type>::value) and
  (::pressio::is_native_container_kokkos<y_type>::value or
    ::pressio::is_expression_acting_on_kokkos<y_type>::value) and
  ::pressio::Traits<A_type>::rank == 2 and
  ::pressio::Traits<x_type>::rank == 1 and
  ::pressio::Traits<y_type>::rank == 1
  >
product(::pressio::nontranspose /*unused*/,
	const scalar_type alpha,
	const A_type & A,
	const x_type & x,
	const scalar_type beta,
	y_type & y)
{
  static_assert
    (::pressio::have_matching_execution_space<A_type, x_type, y_type>::value,
     "operands need to have same execution space" );

  assert( y.extent(0) == A.extent(0) );
  assert( A.extent(1) == x.extent(0) );
  const char ctA = 'N';

  auto & y_n = impl::get_native(y);
  const auto & A_n = impl::get_native(A);
  const auto & x_n = impl::get_native(x);
  ::KokkosBlas::gemv( &ctA, alpha, A_n, x_n, beta, y_n);
}

//-------------------------------
// specialize for op(A) = A^T
//-------------------------------
template < typename A_type, typename x_type, typename scalar_type, typename y_type>
::pressio::mpl::enable_if_t<
  (::pressio::is_native_container_kokkos<A_type>::value or
   ::pressio::is_expression_acting_on_kokkos<A_type>::value) and
  (::pressio::is_native_container_kokkos<x_type>::value or
   ::pressio::is_expression_acting_on_kokkos<x_type>::value) and
  (::pressio::is_native_container_kokkos<y_type>::value or
    ::pressio::is_expression_acting_on_kokkos<y_type>::value) and
  ::pressio::Traits<A_type>::rank == 2 and
  ::pressio::Traits<x_type>::rank == 1 and
  ::pressio::Traits<y_type>::rank == 1
  >
product(::pressio::transpose /*unused*/,
	const scalar_type alpha,
	const A_type & A,
	const x_type & x,
	const scalar_type beta,
	y_type & y)
{
  static_assert
    (::pressio::have_matching_execution_space<A_type,x_type, y_type>::value,
     "operands need to have same execution space" );

  assert( y.extent(0) == A.extent(1) );
  assert( A.extent(0) == x.extent(0) );
  const char ctA = 'T';

  auto & y_n = impl::get_native(y);
  const auto & A_n = impl::get_native(A);
  const auto & x_n = impl::get_native(x);
  ::KokkosBlas::gemv( &ctA, alpha, A_n, x_n, beta, y_n);
}

}}//end namespace pressio::ops
#endif  // OPS_KOKKOS_OPS_LEVEL2_HPP_
