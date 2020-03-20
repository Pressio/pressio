/*
//@HEADER
// ************************************************************************
//
// ops_mvec_prod_mvec.hpp
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

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#ifndef OPS_SRC_OPS_TPETRA_MULTI_VECTOR_PROD_MULTI_VECTOR_HPP_
#define OPS_SRC_OPS_TPETRA_MULTI_VECTOR_PROD_MULTI_VECTOR_HPP_

namespace pressio{ namespace ops{

/*
 * for tpetra:
 *
 * C = beta * C + alpha*op(A)*op(B)
 *
*/

/* -------------------------------------------------------------------
 * specialize for op(A) = A^T and op(B) = B
 *-------------------------------------------------------------------*/
template <
  typename A_type, typename B_type, typename scalar_type, typename C_type
  >
::pressio::mpl::enable_if_t<
  ::pressio::containers::meta::is_multi_vector_wrapper_tpetra<A_type>::value and
  ::pressio::containers::meta::is_multi_vector_wrapper_tpetra<B_type>::value
  /*and ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<C_type>::value*/
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const A_type & A,
	const B_type & B,
	const scalar_type beta,
	::pressio::containers::MatrixSharedMemBase<C_type> & C)
{
  static_assert(containers::meta::are_scalar_compatible<A_type, B_type, C_type>::value,
		"Types are not scalar compatible");

  const auto numVecsA = A.numVectors();
  const auto numVecsB = B.numVectors();
  assert( A.extent(0) == B.extent(0));
  assert(C.extent(0) == numVecsA);
  assert(C.extent(1) == numVecsB);

  // compute dot between every column of A with every col of B
  for (std::size_t i=0; i<(std::size_t)numVecsA; i++)
  {
    // colI is a Teuchos::RCP<Vector<...>>
    const auto colI = A.data()->getVector(i);
    for (std::size_t j=0; j< (std::size_t)numVecsB; j++)
    {
      const auto colJ = B.data()->getVector(j);
      C(i,j) = beta * C(i,j) + alpha * colI->dot(*colJ);
    }
  }
}

template <
  typename C_type, typename A_type, typename B_type, typename scalar_type
  >
::pressio::mpl::enable_if_t<
  ::pressio::containers::meta::is_multi_vector_wrapper_tpetra<A_type>::value and
  ::pressio::containers::meta::is_multi_vector_wrapper_tpetra<B_type>::value and
  ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<C_type>::value,
  C_type
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const A_type & A,
	const B_type & B)
{
  static_assert(containers::meta::are_scalar_compatible<A_type, B_type, C_type>::value,
		"Types are not scalar compatible");
  constexpr auto zero = ::pressio::utils::constants::zero<scalar_type>();

  const auto numVecsA = A.numVectors();
  const auto numVecsB = B.numVectors();
  C_type C(numVecsA, numVecsB);
  product(modeA, modeB, alpha, A, B, zero, C);
  return C;
}


// /***********************************
//  * special case A==B
// **********************************/
template <typename A_type, typename scalar_type, typename C_type>
::pressio::mpl::enable_if_t<
  ::pressio::containers::meta::is_multi_vector_wrapper_tpetra<A_type>::value and
  ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<C_type>::value
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const A_type & A,
	const scalar_type beta,
	C_type & C)
{
  static_assert(containers::meta::are_scalar_compatible<A_type, C_type>::value,
		"Types are not scalar compatible");

  // how many vectors are in A and mvB
  const auto numVecsA = A.numVectors();
  assert(C.extent(0) == numVecsA);
  assert(C.extent(1) == numVecsA);

  scalar_type tmp = ::pressio::utils::constants::zero<scalar_type>();

  // A dot A = A^T*A, which yields a symmetric matrix
  // only need to compute half and fill remaining entries accordingly
  for (std::size_t i=0; i<(std::size_t)numVecsA; i++)
  {
    // colI is a Teuchos::RCP<Vector<...>>
    auto colI = A.data()->getVector(i);
    for (std::size_t j=i; j<(std::size_t)numVecsA; j++)
    {
      auto colJ = A.data()->getVector(j);
      tmp = alpha*colI->dot(*colJ);
      C(i,j) = beta*C(i,j) + tmp;
      C(j,i) = beta*C(j,i) + tmp;
    }
  }
}


template <typename A_type, typename scalar_type, typename C_type>
::pressio::mpl::enable_if_t<
  ::pressio::containers::meta::is_multi_vector_wrapper_tpetra<A_type>::value and
  ::pressio::containers::meta::is_dense_matrix_wrapper_kokkos<C_type>::value
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const A_type & A,
	const scalar_type beta,
	C_type & C)
{
  static_assert(containers::meta::are_scalar_compatible<A_type, C_type>::value,
		"Types are not scalar compatible");

  static_assert(std::is_same<
		typename containers::details::traits<A_type>::device_t,
		typename containers::details::traits<C_type>::device_t
		>::value,
		"Non-matching device types");

  using map_t	    = typename ::pressio::containers::details::traits<A_type>::data_map_t;
  using tpetra_mv_t = typename ::pressio::containers::details::traits<A_type>::wrapped_t;
  const auto indexBase = A.data()->getMap()->getIndexBase();
  const auto comm = A.data()->getMap()->getComm();

  // C should be square matrix
  assert( C.extent(0) == C.extent(1) );
  const auto n = C.extent(0);
  Teuchos::RCP<const map_t> replMap(new map_t(n, indexBase, comm, Tpetra::LocallyReplicated));
  // create multivector that views the Kokkos matrix
  tpetra_mv_t Cmv(replMap, *C.data());

  // do the operation
  Cmv.multiply(Teuchos::ETransp::TRANS, Teuchos::ETransp::NO_TRANS, alpha, *A.data(), *A.data(), beta);
}


template <typename C_type, typename A_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  ::pressio::containers::meta::is_multi_vector_wrapper_tpetra<A_type>::value and
  (::pressio::containers::meta::is_dynamic_dense_matrix_wrapper_eigen<C_type>::value or
   ::pressio::containers::meta::is_dense_matrix_wrapper_kokkos<C_type>::value),
  C_type
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const A_type & A)
{
  static_assert(containers::meta::are_scalar_compatible<A_type, C_type>::value,
		"Types are not scalar compatible");

  constexpr auto zero = ::pressio::utils::constants::zero<scalar_type>();
  C_type C(A.numVectors(), A.numVectors());
  product(modeA, modeB, alpha, A, zero, C);
  return C;
}

}}//end namespace pressio::ops
#endif
#endif
