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
C is an Eigen dense matrix with col major
*-------------------------------------------------------------------*/
#ifdef PRESSIO_ENABLE_TPL_EIGEN
template <
  class A_type, class B_type, class C_type,
  class alpha_t, class beta_t
  >
std::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<B_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra<A_type>::value
  && ::pressio::is_multi_vector_tpetra<B_type>::value
  && ::pressio::is_dense_col_major_matrix_eigen<C_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, B_type, C_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t,  typename ::pressio::Traits<A_type>::scalar_type>::value
  >
product(::pressio::transpose /*unused*/,
	::pressio::nontranspose /*unused*/,
	const alpha_t & alpha,
	const A_type & A,
	const B_type & B,
	const beta_t & beta,
	C_type & C)
{

  assert( (std::size_t)::pressio::ops::extent(A,0) == (std::size_t)::pressio::ops::extent(B,0));
  assert( (std::size_t)::pressio::ops::extent(C,0) == (std::size_t) A.getNumVectors() );
  assert( (std::size_t)::pressio::ops::extent(C,1) == (std::size_t) B.getNumVectors() );

  // for col-major, we should have inner stride == 1 and outerstride = num of rows
  // we should actually check at compile time somehow if C can be viewed in a kokkos view
  assert(C.innerStride() == 1);
  assert(C.outerStride() == C.rows());

  using C_sc_t = typename ::pressio::Traits<C_type>::scalar_type;
  using kokkos_view_t = Kokkos::View<C_sc_t**, Kokkos::LayoutLeft, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
  kokkos_view_t Cview(C.data(), C.rows(), C.cols());

  const C_sc_t alpha_(alpha);
  const C_sc_t beta_(beta);

  using map_t = typename A_type::map_type;
  const auto indexBase = A.getMap()->getIndexBase();
  const auto comm = A.getMap()->getComm();
  Teuchos::RCP<const map_t> replMap(new map_t(Cview.extent(0), indexBase,
					      comm, Tpetra::LocallyReplicated));
  A_type Cmv(replMap, Cview);
  Cmv.multiply(Teuchos::ETransp::TRANS, Teuchos::ETransp::NO_TRANS, alpha_, A, B, beta_);


  /* this is the old code used when using the Eigen matrix directly
     but this was inefficient.

     // compute dot between every column of A with every col of B
     for (std::size_t i=0; i<(std::size_t)numVecsA; i++){
       // colI is a Teuchos::RCP<Vector<...>>
       const auto colI = A.getVector(i);
       for (std::size_t j=0; j<(std::size_t)numVecsB; j++)
       {
	 const auto colJ = B.getVector(j);
	 if (beta == static_cast<beta_t>(0)) {
	   C(i,j) = alpha * colI->dot(*colJ);
	 } else {
	   C(i,j) = beta * C(i,j) + alpha * colI->dot(*colJ);
	 }
       }
     }
  */
}

/* -------------------------------------------------------------------
C = beta * C + alpha * A^T * B

A = tpetra multivector
B = tpetra multivector
C is an Eigen dense matrix returned by the function
*-------------------------------------------------------------------*/
template <
  class C_type, class A_type, class B_type,
  class alpha_t
  >
std::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<B_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra<A_type>::value
  && ::pressio::is_multi_vector_tpetra<B_type>::value
  && ::pressio::is_dynamic_dense_matrix_eigen<C_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, B_type, C_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value,
  C_type
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const alpha_t & alpha,
	const A_type & A,
	const B_type & B)
{

  using sc_t = typename ::pressio::Traits<C_type>::scalar_type;
  const auto numVecsA = A.getNumVectors();
  const auto numVecsB = B.getNumVectors();
  C_type C(numVecsA, numVecsB);
  product(modeA, modeB, alpha, A, B, static_cast<sc_t>(0), C);
  return C;
}

// /***********************************
//  * special case A==B
// **********************************/
template <
  class A_type, class C_type,
  class alpha_t, class beta_t
  >
std::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra<A_type>::value
  && ::pressio::is_dense_matrix_eigen<C_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, C_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t,  typename ::pressio::Traits<A_type>::scalar_type>::value
  >
product(::pressio::transpose /*unused*/,
	::pressio::nontranspose /*unused*/,
	const alpha_t & alpha,
	const A_type & A,
	const beta_t & beta,
	C_type & C)
{

  // how many vectors are in A and mvB
  const auto numVecsA = A.getNumVectors();
  assert((std::size_t)C.rows() == (std::size_t)numVecsA);
  assert((std::size_t)C.cols() == (std::size_t)numVecsA);
  using sc_t = typename ::pressio::Traits<A_type>::scalar_type;
  const auto zero = static_cast<sc_t>(0);
  const sc_t alpha_(alpha);
  const sc_t beta_(beta);
  const auto has_beta = beta_ != zero;
  sc_t tmp = zero;

  // A dot A = A^T*A, which yields a symmetric matrix
  // only need to compute half and fill remaining entries accordingly
  for (std::size_t i=0; i<(std::size_t)numVecsA; i++)
  {
    // colI is a Teuchos::RCP<Vector<...>>
    auto colI = A.getVector(i);
    for (std::size_t j=i; j<(std::size_t)numVecsA; j++)
    {
      auto colJ = A.getVector(j);
      tmp = alpha_*colI->dot(*colJ);
      C(i,j) = has_beta ? beta_*C(i,j) + tmp : tmp;
      if(j!=i){
        C(j,i) = has_beta ? beta_*C(j,i) + tmp : tmp;
      }
    }
  }
}

template <class C_type, class A_type, class alpha_t>
std::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra<A_type>::value
  && ::pressio::is_dynamic_dense_matrix_eigen<C_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, C_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value,
  C_type
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const alpha_t & alpha,
	const A_type & A)
{

  using sc_t = typename ::pressio::Traits<C_type>::scalar_type;
  C_type C(A.getNumVectors(), A.getNumVectors());
  product(modeA, modeB, alpha, A, static_cast<sc_t>(0), C);
  return C;
}
#endif


/* -------------------------------------------------------------------
C = beta * C + alpha * A^T * A

A = tpetra multivector
C is a Kokkos dense matrix
*-------------------------------------------------------------------*/
template <class A_type, class C_type, class alpha_t, class beta_t>
std::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra<A_type>::value
  && ::pressio::is_dense_matrix_kokkos<C_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, C_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  >
product(::pressio::transpose /*unused*/,
	::pressio::nontranspose /*unused*/,
	const alpha_t & alpha,
	const A_type & A,
	const beta_t & beta,
	C_type & C)
{

  static_assert
    (std::is_same< typename C_type::array_layout, Kokkos::LayoutLeft>::value,
     "The kokkos matrix must be layout left");

  assert( (std::size_t)::pressio::ops::extent(C, 0) == (std::size_t) A.getNumVectors() );
  assert( (std::size_t)::pressio::ops::extent(C, 1) == (std::size_t) A.getNumVectors() );

  using map_t      = typename A_type::map_type;
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
  using sc_t = typename ::pressio::Traits<A_type>::scalar_type;
  const sc_t alpha_(alpha);
  const sc_t beta_(beta);
  Cmv.multiply(Teuchos::ETransp::TRANS, Teuchos::ETransp::NO_TRANS, alpha_, A, A, beta_);
}

/* -------------------------------------------------------------------
C = alpha * A^T * A

A = tpetra multivector
C is a Kokkos dense matrix
*-------------------------------------------------------------------*/
template <class C_type, class A_type, class alpha_t>
std::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra<A_type>::value
  && ::pressio::is_dynamic_dense_matrix_kokkos<C_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, C_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value,
  C_type
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const alpha_t & alpha,
	const A_type & A)
{

  using sc_t = typename ::pressio::Traits<C_type>::scalar_type;
  C_type C("opsProdC", A.getNumVectors(), A.getNumVectors());
  product(modeA, modeB, alpha, A, static_cast<sc_t>(0), C);
  return C;
}

/* -------------------------------------------------------------------
C = beta * C + alpha * A^T * B

A = tpetra multivector
B = tpetra multivector
C is a Kokkos dense matrix
*-------------------------------------------------------------------*/
template <class A_type, class B_type, class C_type, class alpha_t, class beta_t>
std::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<B_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra<A_type>::value
  && ::pressio::is_multi_vector_tpetra<B_type>::value
  && ::pressio::is_dense_matrix_kokkos<C_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, B_type, C_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t,  typename ::pressio::Traits<A_type>::scalar_type>::value
  >
product(::pressio::transpose /*unused*/,
	::pressio::nontranspose /*unused*/,
	const alpha_t & alpha,
	const A_type & A,
	const B_type & B,
	const beta_t & beta,
	C_type & C)
{
  static_assert
    (std::is_same< typename C_type::array_layout, Kokkos::LayoutLeft>::value,
     "The kokkos matrix must be layout left");

  assert( (std::size_t)::pressio::ops::extent(A, 0) == (std::size_t)::pressio::ops::extent(B, 0));
  assert( (std::size_t)::pressio::ops::extent(C, 0) == (std::size_t) A.getNumVectors() );
  assert( (std::size_t)::pressio::ops::extent(C, 1) == (std::size_t) B.getNumVectors() );

  using map_t = typename A_type::map_type;
  const auto indexBase = A.getMap()->getIndexBase();
  const auto comm = A.getMap()->getComm();

  const auto n = C.extent(0);
  Teuchos::RCP<const map_t> replMap(new map_t(n, indexBase, comm, Tpetra::LocallyReplicated));
  // create multivector that views the Kokkos matrix
  A_type Cmv(replMap, C);

  // do the operation
  using sc_t = typename ::pressio::Traits<A_type>::scalar_type;
  const sc_t alpha_(alpha);
  const sc_t beta_(beta);
  Cmv.multiply(Teuchos::ETransp::TRANS, Teuchos::ETransp::NO_TRANS, alpha_, A, B, beta_);
}

}}//end namespace pressio::ops





// //-------------------------------------------
// // C = beta * C + alpha*A*B
// // specialize for when A is a diagonal expression
// //-------------------------------------------
// template <typename T, typename B_type, typename scalar_type, typename C_type>
// std::enable_if_t<
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
//   assert( C.extent(0) == A.extent(0) );
//   assert( C.extent(1) == B.extent(1) );
//   assert( A.extent(1) == B.extent(0) );

//   auto & Ctp = *C.data();
//   auto & Atp = *(A.pressioObj()->data());
//   auto & Btp = *B.data();
//   Ctp.elementWiseMultiply(alpha, Atp, Btp, beta);
// }
#endif  // OPS_TPETRA_OPS_LEVEL3_HPP_
