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

#ifndef OPS_TPETRA_OPS_LEVEL3_HPP_
#define OPS_TPETRA_OPS_LEVEL3_HPP_

namespace pressio{ namespace ops{

/* -------------------------------------------------------------------
C = beta * C + alpha * A^T * B

A = tpetra multivector
B = tpetra multivector
C is an Eigen dense matrix
*-------------------------------------------------------------------*/
#ifdef PRESSIO_ENABLE_TPL_EIGEN
template <
  typename A_type, typename B_type, typename scalar_type, typename C_type
  >
::pressio::mpl::enable_if_t<
  ::pressio::is_multi_vector_tpetra<A_type>::value and
  ::pressio::is_multi_vector_tpetra<B_type>::value and
  ::pressio::is_dense_matrix_eigen<C_type>::value
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const A_type & A,
	const B_type & B,
	const scalar_type beta,
	C_type & C)
{
  static_assert
    (are_scalar_compatible<A_type, B_type, C_type>::value,
     "Types are not scalar compatible");

  const auto numVecsA = A.getNumVectors();
  const auto numVecsB = B.getNumVectors();
  assert( (std::size_t)::pressio::ops::extent(A,0) == (std::size_t)::pressio::ops::extent(B,0));
  assert( (std::size_t)::pressio::ops::extent(C,0) == (std::size_t)numVecsA );
  assert( (std::size_t)::pressio::ops::extent(C,1) == (std::size_t)numVecsB );

  // compute dot between every column of A with every col of B
  for (std::size_t i=0; i<(std::size_t)numVecsA; i++)
  {
    // colI is a Teuchos::RCP<Vector<...>>
    const auto colI = A.getVector(i);
    for (std::size_t j=0; j<(std::size_t)numVecsB; j++)
    {
      const auto colJ = B.getVector(j);
      C(i,j) = beta * C(i,j) + alpha * colI->dot(*colJ);
    }
  }
}

/* -------------------------------------------------------------------
C = beta * C + alpha * A^T * B

A = tpetra multivector
B = tpetra multivector
C is an Eigen dense matrix returned by the function
*-------------------------------------------------------------------*/
template <
  typename C_type, typename A_type, typename B_type, typename scalar_type
  >
::pressio::mpl::enable_if_t<
  ::pressio::is_multi_vector_tpetra<A_type>::value and
  ::pressio::is_multi_vector_tpetra<B_type>::value and
  ::pressio::is_dense_matrix_eigen<C_type>::value,
  C_type
  >
product(::pressio::transpose modeA,
 ::pressio::nontranspose modeB,
 const scalar_type alpha,
 const A_type & A,
 const B_type & B)
{
  static_assert
    (are_scalar_compatible<A_type, B_type, C_type>::value,
     "Types are not scalar compatible");
  constexpr auto zero = ::pressio::utils::Constants<scalar_type>::zero();

  const auto numVecsA = A.getNumVectors();
  const auto numVecsB = B.getNumVectors();
  C_type C(numVecsA, numVecsB);
  product(modeA, modeB, alpha, A, B, zero, C);
  return C;
}

// /***********************************
//  * special case A==B
// **********************************/
template <typename A_type, typename scalar_type, typename C_type>
::pressio::mpl::enable_if_t<
  ::pressio::is_multi_vector_tpetra<A_type>::value and
  ::pressio::is_dense_matrix_eigen<C_type>::value
  >
product(::pressio::transpose modeA,
 ::pressio::nontranspose modeB,
 const scalar_type alpha,
 const A_type & A,
 const scalar_type beta,
 C_type & C)
{
  static_assert(are_scalar_compatible<A_type, C_type>::value,
   "Types are not scalar compatible");

  // how many vectors are in A and mvB
  const auto numVecsA = A.getNumVectors();
  assert((std::size_t)C.rows() == (std::size_t)numVecsA);
  assert((std::size_t)C.cols() == (std::size_t)numVecsA);
  scalar_type tmp = ::pressio::utils::Constants<scalar_type>::zero();

  // A dot A = A^T*A, which yields a symmetric matrix
  // only need to compute half and fill remaining entries accordingly
  for (std::size_t i=0; i<(std::size_t)numVecsA; i++)
  {
    // colI is a Teuchos::RCP<Vector<...>>
    auto colI = A.getVector(i);
    for (std::size_t j=i; j<(std::size_t)numVecsA; j++)
    {
      auto colJ = A.getVector(j);
      tmp = alpha*colI->dot(*colJ);
      C(i,j) = beta*C(i,j) + tmp;
      if(j!=i){
        C(j,i) = beta*C(j,i) + tmp;
      }
    }
  }
}

template <typename C_type, typename A_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  ::pressio::is_multi_vector_tpetra<A_type>::value and
  ::pressio::is_dynamic_dense_matrix_eigen<C_type>::value,
  C_type
  >
product(::pressio::transpose modeA,
 ::pressio::nontranspose modeB,
 const scalar_type alpha,
 const A_type & A)
{
  static_assert(are_scalar_compatible<A_type, C_type>::value,
   "Types are not scalar compatible");

  constexpr auto zero = ::pressio::utils::Constants<scalar_type>::zero();
  C_type C(A.getNumVectors(), A.getNumVectors());
  product(modeA, modeB, alpha, A, zero, C);
  return C;
}
#endif


/* -------------------------------------------------------------------
C = beta * C + alpha * A^T * A

A = tpetra multivector
C is a Kokkos dense matrix
*-------------------------------------------------------------------*/
template <typename A_type, typename scalar_type, typename C_type>
::pressio::mpl::enable_if_t<
  ::pressio::is_multi_vector_tpetra<A_type>::value and
  ::pressio::is_dense_matrix_kokkos<C_type>::value
  >
product(::pressio::transpose modeA,
 ::pressio::nontranspose modeB,
 const scalar_type alpha,
 const A_type & A,
 const scalar_type beta,
 C_type & C)
{
  static_assert
    (are_scalar_compatible<A_type, C_type>::value,
     "Types are not scalar compatible");

  static_assert
    (std::is_same<
     typename ::pressio::Traits<A_type>::device_type,
     typename ::pressio::Traits<C_type>::device_type
     >::value,
     "Non-matching device types");

  static_assert
    (std::is_same< typename C_type::array_layout, Kokkos::LayoutLeft>::value,
     "The kokkos matrix must be layout left");

  using map_t      = typename ::pressio::Traits<A_type>::data_map_type;
  // using tpetra_mv_t = typename ::pressio::Traits<A_type>::wrapped_type;
  const auto indexBase = A.getMap()->getIndexBase();
  const auto comm = A.getMap()->getComm();

  // C should be square matrix
  assert(C.extent(0) == C.extent(1));
  const auto n = C.extent(0);
  Teuchos::RCP<const map_t> replMap(new map_t(n, indexBase, comm, Tpetra::LocallyReplicated));
  // create multivector that views the Kokkos matrix
  A_type Cmv(replMap, C);

  // do the operation
  Cmv.multiply(Teuchos::ETransp::TRANS, Teuchos::ETransp::NO_TRANS, alpha, A, A, beta);
}

/* -------------------------------------------------------------------
C = alpha * A^T * A

A = tpetra multivector
C is a Kokkos dense matrix
*-------------------------------------------------------------------*/
template <class C_type, typename A_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  ::pressio::is_multi_vector_tpetra<A_type>::value and
  ::pressio::is_dense_matrix_kokkos<C_type>::value,
  C_type
  >
product(::pressio::transpose modeA,
 ::pressio::nontranspose modeB,
 const scalar_type alpha,
 const A_type & A)
{
  static_assert(are_scalar_compatible<A_type, C_type>::value,
    "Types are not scalar compatible");

  constexpr auto zero = ::pressio::utils::Constants<scalar_type>::zero();
  C_type C("opsProdC", A.getNumVectors(), A.getNumVectors());
  product(modeA, modeB, alpha, A, zero, C);
  return C;
}


/* -------------------------------------------------------------------
C = beta * C + alpha * A^T * B

A = tpetra multivector
B = tpetra multivector
C is a Kokkos dense matrix
*-------------------------------------------------------------------*/
template <typename A_type, typename scalar_type, typename C_type>
::pressio::mpl::enable_if_t<
  ::pressio::is_multi_vector_tpetra<A_type>::value and
  ::pressio::is_dense_matrix_kokkos<C_type>::value
  >
product(::pressio::transpose modeA,
 ::pressio::nontranspose modeB,
 const scalar_type alpha,
 const A_type & A,
 const A_type & B,
 const scalar_type beta,
 C_type & C)
{
  static_assert
    (are_scalar_compatible<A_type, C_type>::value,
     "Types are not scalar compatible");

  static_assert
    (std::is_same<
     typename ::pressio::Traits<A_type>::device_type,
     typename ::pressio::Traits<C_type>::device_type
     >::value,
     "Non-matching device types");

  static_assert
    (std::is_same< typename C_type::array_layout, Kokkos::LayoutLeft>::value,
     "The kokkos matrix must be layout left");

  using map_t      = typename ::pressio::Traits<A_type>::data_map_type;
  const auto indexBase = A.getMap()->getIndexBase();
  const auto comm = A.getMap()->getComm();

  const auto n = C.extent(0);
  Teuchos::RCP<const map_t> replMap(new map_t(n, indexBase, comm, Tpetra::LocallyReplicated));
  // create multivector that views the Kokkos matrix
  A_type Cmv(replMap, C);

  // do the operation
  Cmv.multiply(Teuchos::ETransp::TRANS, Teuchos::ETransp::NO_TRANS, alpha, A, B, beta);
}

}}//end namespace pressio::ops




// //-------------------------------------------
// // C = beta * C + alpha*A*B
// // specialize for when A is a diagonal expression
// //-------------------------------------------
// template <typename T, typename B_type, typename scalar_type, typename C_type>
// ::pressio::mpl::enable_if_t<
//   ::pressio::containers::predicates::is_multi_vector_wrapper_tpetra<B_type>::value and
//   ::pressio::containers::predicates::is_multi_vector_wrapper_tpetra<C_type>::value and
//   ::pressio::containers::predicates::is_vector_wrapper_tpetra<T>::value
//   >
// product(::pressio::nontranspose modeA,
// 	::pressio::nontranspose modeB,
// 	const scalar_type alpha,
// 	const pressio::containers::expressions::AsDiagonalMatrixExpr<T> & A,
// 	const B_type & B,
// 	const scalar_type beta,
// 	C_type & C)
// {
//   static_assert
//     (containers::predicates::are_scalar_compatible<T, B_type, C_type>::value,
//      "Types are not scalar compatible");

//   assert( C.extent(0) == A.extent(0) );
//   assert( C.extent(1) == B.extent(1) );
//   assert( A.extent(1) == B.extent(0) );

//   auto & Ctp = *C.data();
//   auto & Atp = *(A.pressioObj()->data());
//   auto & Btp = *B.data();
//   Ctp.elementWiseMultiply(alpha, Atp, Btp, beta);
// }
#endif  // OPS_TPETRA_OPS_LEVEL3_HPP_
