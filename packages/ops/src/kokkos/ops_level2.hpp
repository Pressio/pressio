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
  ::pressio::ops::constraints::rank2_container_kokkos_with_native_data_access<A_type>::value and
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<x_type>::value
  >
product(::pressio::nontranspose mode,
	const scalar_type alpha,
	const A_type & A,
	const x_type & x,
	const scalar_type beta,
	y_type & y)
{
  /* make sure we don't pass const objects to be modified.
     In kokkos it is legal to modify const views, not for pressio wrappers. */
  static_assert
    (!std::is_const<y_type>::value,
     "ops:product: cannot modify a const-qualified wrapper of a Kokkos view");
  static_assert
    (containers::predicates::are_scalar_compatible<A_type, x_type, y_type>::value,
     "Types are not scalar compatible");
  static_assert
    (::pressio::containers::predicates::have_matching_execution_space<
     A_type, x_type, y_type>::value,
     "operands need to have same execution space" );

  assert( y.extent(0) == A.extent(0) );
  assert( A.extent(1) == x.extent(0) );
  const char ctA = 'N';
  ::KokkosBlas::gemv( &ctA, alpha, *A.data(), *x.data(), beta, *y.data() );
}

//-------------------------------
// specialize for op(A) = A^T
//-------------------------------
template < typename A_type, typename x_type, typename scalar_type, typename y_type>
::pressio::mpl::enable_if_t<
  ::pressio::ops::constraints::rank2_container_kokkos_with_native_data_access<A_type>::value and
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<x_type>::value and
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<y_type>::value
  >
product(::pressio::transpose mode,
	const scalar_type alpha,
	const A_type & A,
	const x_type & x,
	const scalar_type beta,
	y_type & y)
{
  /* make sure we don't pass const objects to be modified.
     In kokkos it is legal to modify const views, not for pressio wrappers.
   */
  static_assert
    (!std::is_const<y_type>::value,
     "ops:product: cannot modify a const-qualified wrapper of a Kokkos view");

  static_assert
    (containers::predicates::are_scalar_compatible<A_type, x_type, y_type>::value,
     "Types are not scalar compatible");
  static_assert
    (::pressio::containers::predicates::have_matching_execution_space<
     A_type,x_type, y_type>::value,
     "operands need to have same execution space" );

  assert( y.extent(0) == A.extent(1) );
  assert( A.extent(0) == x.extent(0) );
  const char ctA = 'T';
  ::KokkosBlas::gemv( &ctA, alpha, *A.data(), *x.data(), beta, *y.data() );
}

}}//end namespace pressio::ops
#endif  // OPS_KOKKOS_OPS_LEVEL2_HPP_
