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
 * for tpetra:
 *
 * C = beta * C + alpha*op(A)*op(B)
 *
*/

//-------------------------------------------
// specialize for op(A) = A^T and op(B) = B
//-------------------------------------------
template <typename A_type, typename B_type, typename scalar_type, typename C_type>
::pressio::mpl::enable_if_t<
  ::pressio::ops::constraints::rank2_container_kokkos_with_native_data_access<A_type>::value and
  ::pressio::ops::constraints::rank2_container_kokkos_with_native_data_access<B_type>::value and
  ::pressio::ops::constraints::rank2_container_kokkos_with_native_data_access<C_type>::value
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const A_type & A,
	const B_type & B,
	const scalar_type beta,
	C_type & C)
{
  /* make sure we don't pass const objects to be modified.
     In kokkos it is legal to modify const views, not for pressio wrappers. */
  static_assert
    (!std::is_const<C_type>::value,
     "ops:product: cannot modify a const-qualified wrapper of a Kokkos view");
  static_assert
    (containers::predicates::are_scalar_compatible<A_type, B_type, C_type>::value,
     "Types are not scalar compatible");
  static_assert
    (::pressio::containers::predicates::have_matching_execution_space<
     A_type, B_type, C_type>::value,
     "operands need to have same execution space" );

  const char ctA = 'T';
  const char ctB = 'N';
  ::KokkosBlas::gemm(&ctA, &ctB, alpha, *A.data(), *B.data(), beta, *C.data());
}

/*----------------------------------------
* special case A==B and op(A)=A^T, op(B)=B
------------------------------------------*/
template <typename A_type, typename scalar_type, typename C_type>
::pressio::mpl::enable_if_t<
  ::pressio::ops::constraints::rank2_container_kokkos_with_native_data_access<A_type>::value and
  ::pressio::ops::constraints::rank2_container_kokkos_with_native_data_access<C_type>::value
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
  ::pressio::ops::constraints::rank2_container_kokkos_with_native_data_access<A_type>::value and
  ::pressio::ops::constraints::rank2_container_kokkos_with_native_data_access<C_type>::value,
  C_type
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const A_type & A)
{
  constexpr auto zero = ::pressio::utils::constants<scalar_type>::zero();
  C_type C(A.numVectors(), A.numVectors());
  product(modeA, modeB, alpha, A, A, zero, C);
  return C;
}

}}//end namespace pressio::ops
#endif  // OPS_KOKKOS_OPS_LEVEL3_HPP_
