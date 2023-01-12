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

#ifndef OPS_TPETRA_BLOCK_OPS_LEVEL3_HPP_
#define OPS_TPETRA_BLOCK_OPS_LEVEL3_HPP_

namespace pressio{ namespace ops{

/*
 * C = beta * C + alpha*op(A)*op(B)
*/

#ifdef PRESSIO_ENABLE_TPL_EIGEN
/* -------------------------------------------------------------------
C = beta * C + alpha * A^T * B

A = tpetra block multivector
B = tpetra block multivector
C is an Eigen dense matrix
*-------------------------------------------------------------------*/
template <
  class A_type, class B_type, class C_type, class alpha_t, class beta_t
  >
::pressio::mpl::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<B_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra_block<A_type>::value
  && ::pressio::is_multi_vector_tpetra_block<B_type>::value
  && ::pressio::is_dense_col_major_matrix_eigen<C_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, B_type, C_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t,  typename ::pressio::Traits<A_type>::scalar_type>::value
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const alpha_t & alpha,
	const A_type & A,
	const B_type & B,
	const beta_t & beta,
	C_type & C)
{
  auto A_tp = A.getMultiVectorView();
  auto B_tp = B.getMultiVectorView();
  ::pressio::ops::product(modeA, modeB, alpha, A_tp, B_tp, beta, C);
}

/* -------------------------------------------------------------------
C = beta * C + alpha * A^T * B

A = tpetra block multivector
B = tpetra block multivector
C is an Eigen dense matrix returned by the function
*-------------------------------------------------------------------*/
template <
  class C_type, class A_type, class B_type, class alpha_t
  >
::pressio::mpl::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<B_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra_block<A_type>::value
  && ::pressio::is_multi_vector_tpetra_block<B_type>::value
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

  auto A_tp = A.getMultiVectorView();
  auto B_tp = B.getMultiVectorView();
  return product<C_type>(modeA, modeB, alpha, A_tp, B_tp);
}

// /***********************************
//  * special case A==B
// **********************************/
template <
  class A_type, class C_type, class alpha_t, class beta_t
  >
::pressio::mpl::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra_block<A_type>::value
  && ::pressio::is_dense_matrix_eigen<C_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, C_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const alpha_t & alpha,
	const A_type & A,
	const beta_t & beta,
	C_type & C)
{
  auto A_tp = A.getMultiVectorView();
  product(modeA, modeB, alpha, A_tp, beta, C);
}

template <
  class C_type, class A_type, class alpha_t
  >
::pressio::mpl::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra_block<A_type>::value
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
  auto A_tp = A.getMultiVectorView();
  return product<C_type>(modeA, modeB, alpha, A_tp);
}
#endif


/* -------------------------------------------------------------------
C = beta * C + alpha * A^T * A

A = tpetra block multivector
C is a Kokkos dense matrix
*-------------------------------------------------------------------*/
template <
  class A_type, class C_type, class alpha_t, class beta_t
  >
::pressio::mpl::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra_block<A_type>::value
  && ::pressio::is_dense_matrix_kokkos<C_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, C_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t,  typename ::pressio::Traits<A_type>::scalar_type>::value
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const alpha_t & alpha,
	const A_type & A,
	const beta_t & beta,
	C_type & C)
{
  auto A_tp = A.getMultiVectorView();
  product(modeA, modeB, alpha, A_tp, beta, C);
}

/* -------------------------------------------------------------------
   C = alpha * A^T * A

   A = tpetra block multivector
   C is a Kokkos dense matrix
   *-------------------------------------------------------------------*/
template <class C_type, class A_type, class alpha_t>
::pressio::mpl::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra_block<A_type>::value
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
  auto A_tp = A.getMultiVectorView();
  return product<C_type>(modeA, modeB, alpha, A_tp);
}

/* -------------------------------------------------------------------
C = beta * C + alpha * A^T * B

A = tpetra block multivector
B = tpetra block multivector
C is a Kokkos dense matrix
*-------------------------------------------------------------------*/
template <class A_type, class C_type, class alpha_t, class beta_t>
::pressio::mpl::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra_block<A_type>::value
  && ::pressio::is_dense_matrix_kokkos<C_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, C_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t,  typename ::pressio::Traits<A_type>::scalar_type>::value
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const alpha_t & alpha,
	const A_type & A,
	const A_type & B,
	const beta_t & beta,
	C_type & C)
{
  auto A_tp = A.getMultiVectorView();
  auto B_tp = B.getMultiVectorView();
  product(modeA, modeB, alpha, A_tp, B_tp, beta, C);
}

}}//end namespace pressio::ops









// /* -------------------------------------------------------------------
//  * specialize for op(A) = A^T and op(B) = B
//  *-------------------------------------------------------------------*/
// template <
//   class A_type, class B_type, class scalar_type, class C_type
//   >
// ::pressio::mpl::enable_if_t<
//   ::pressio::containers::predicates::is_multi_vector_wrapper_tpetra_block<A_type>::value and
//   ::pressio::containers::predicates::is_multi_vector_wrapper_tpetra_block<B_type>::value and
//   ::pressio::ops::constraints::sharedmem_host_subscriptable_rank2_container<C_type>::value
//   >
// product(::pressio::transpose modeA,
// 	::pressio::nontranspose modeB,
// 	const scalar_type alpha,
// 	const A_type & A,
// 	const B_type & B,
// 	const scalar_type beta,
// 	C_type & C)
// {
//   // get a tpetra multivector that views the data
//   const auto Amvv = A.data()->getMultiVectorView();
//   const auto Bmvv = B.data()->getMultiVectorView();
//   const auto numVecsA = A.numVectors();
//   const auto numVecsB = B.numVectors();
//   assert((std::size_t)A.extent(0) == (std::size_t)B.extent(0));
//   assert((std::size_t)C.extent(0) == (std::size_t)numVecsA);
//   assert((std::size_t)C.extent(1) == (std::size_t)numVecsB);

//   for (std::size_t i=0; i<(std::size_t)numVecsA; i++)
//   {
//     // colI is a Teuchos::RCP<Vector<...>>
//     auto colI = Amvv.getVector(i);
//     for (std::size_t j=0; j<(std::size_t)numVecsB; j++)
//     {
//       auto colJ = Bmvv.getVector(j);
//       C(i,j) = beta*C(i,j) + alpha*colI->dot(*colJ);
//     }
//   }
// }

// template <
//   class C_type, class A_type, class B_type, class scalar_type
//   >
// ::pressio::mpl::enable_if_t<
//   ::pressio::containers::predicates::is_multi_vector_wrapper_tpetra_block<A_type>::value and
//   ::pressio::containers::predicates::is_multi_vector_wrapper_tpetra_block<B_type>::value and
//   ::pressio::ops::constraints::sharedmem_host_subscriptable_rank2_container<C_type>::value and
//   !::pressio::ops::constraints::rank2_container_kokkos_with_native_data_access<C_type>::value,
//   C_type
//   >
// product(::pressio::transpose modeA,
// 	::pressio::nontranspose modeB,
// 	const scalar_type alpha,
// 	const A_type & A,
// 	const B_type & B)
// {
//   constexpr auto zero = ::pressio::utils::Constants<scalar_type>::zero();

//   const auto numVecsA = A.numVectors();
//   const auto numVecsB = B.numVectors();
//   C_type C(numVecsA, numVecsB);
//   product(modeA, modeB, alpha, A, B, zero, C);
//   return C;
// }

// // /***********************************
// //  * special case A==B
// // **********************************/
// template <class A_type, class scalar_type, class C_type>
// ::pressio::mpl::enable_if_t<
//   ::pressio::containers::predicates::is_multi_vector_wrapper_tpetra_block<A_type>::value and
//   ::pressio::ops::constraints::sharedmem_host_subscriptable_rank2_container<C_type>::value and
//   !::pressio::containers::predicates::is_dense_matrix_wrapper_kokkos<C_type>::value
//   >
// product(::pressio::transpose modeA,
// 	::pressio::nontranspose modeB,
// 	const scalar_type alpha,
// 	const A_type & A,
// 	const scalar_type beta,
// 	C_type & C)
// {
//   // get a tpetra multivector that views the data
//   const auto Amvv = A.data()->getMultiVectorView();
//   const auto numVecsA = A.numVectors();
//   assert(C.extent(0) == numVecsA);
//   assert(C.extent(1) == numVecsA);
//   scalar_type tmp = ::pressio::utils::Constants<scalar_type>::zero();

//   // using ord_t = typename ::pressio::containers::details::traits<A_type>::global_ordinal_t;

//   // A dot A = A^T*A, which yields a symmetric matrix
//   // only need to compute half and fill remaining entries accordingly
//   for (std::size_t i=0; i<(std::size_t)numVecsA; i++)
//   {
//     // colI is a Teuchos::RCP<Vector<...>>
//     auto colI = Amvv.getVector(i);
//     for (std::size_t j=i; j<(std::size_t)numVecsA; j++)
//     {
//       auto colJ = Amvv.getVector(j);
//       tmp = alpha*colI->dot(*colJ);
//       C(i,j) = beta*C(i,j) + tmp;
//       C(j,i) = beta*C(j,i) + tmp;
//     }
//   }
// }

// template <class A_type, class scalar_type, class C_type>
// ::pressio::mpl::enable_if_t<
//   ::pressio::containers::predicates::is_multi_vector_wrapper_tpetra_block<A_type>::value and
//   ::pressio::ops::constraints::rank2_container_kokkos_with_native_data_access<C_type>::value
//   >
// product(::pressio::transpose modeA,
// 	::pressio::nontranspose modeB,
// 	const scalar_type alpha,
// 	const A_type & A,
// 	const scalar_type beta,
// 	C_type & C)
// {
//   // check traits of the block mv
//   using tpetra_mvb_t = typename ::pressio::containers::details::traits<A_type>::wrapped_t;

//   // from the mvb type we can get the underlying regular tpetra::mv and map
//   using tpetra_mv_t  = typename tpetra_mvb_t::mv_type;
//   using map_t	     = typename tpetra_mv_t::map_type;

//   // tpetra multivector that views the tpetra block data
//   const auto mvView = A.data()->getMultiVectorView();

//   const auto indexBase = mvView.getMap()->getIndexBase();
//   const auto comm = mvView.getMap()->getComm();
//   // C should be symmetric
//   assert( (std::size_t)C.extent(0) == (std::size_t)C.extent(1) );
//   const auto n = C.extent(0);
//   Teuchos::RCP<const map_t> replMap(new map_t(n, indexBase, comm, Tpetra::LocallyReplicated));
//   // create multivector that views the Kokkos matrix result
//   tpetra_mv_t Cmv(replMap, *C.data());

//   // do the operation C = A^T A
//   Cmv.multiply(Teuchos::ETransp::TRANS, Teuchos::ETransp::NO_TRANS,
// 	       alpha, mvView, mvView, beta);
// }

// template <class C_type, class A_type, class scalar_type>
// ::pressio::mpl::enable_if_t<
//   ::pressio::containers::predicates::is_multi_vector_wrapper_tpetra_block<A_type>::value and
//   (::pressio::ops::constraints::sharedmem_host_subscriptable_rank2_container<C_type>::value or
//    ::pressio::ops::constraints::rank2_container_kokkos_with_native_data_access<C_type>::value),
//   C_type
//   >
// product(::pressio::transpose modeA,
// 	::pressio::nontranspose modeB,
// 	const scalar_type alpha,
// 	const A_type & A)
// {
//   constexpr auto zero = ::pressio::utils::Constants<scalar_type>::zero();
//   C_type C(A.numVectors(), A.numVectors());
//   product(modeA, modeB, alpha, A, zero, C);
//   return C;
// }

// //-------------------------------------------
// // C = beta * C + alpha*A*B
// // specialize for when A is a diagonal expression
// //-------------------------------------------
// template <class T, class B_type, class scalar_type, class C_type>
// ::pressio::mpl::enable_if_t<
//   ::pressio::containers::predicates::is_multi_vector_wrapper_tpetra_block<B_type>::value and
//   ::pressio::containers::predicates::is_multi_vector_wrapper_tpetra_block<C_type>::value and
//   ::pressio::containers::predicates::is_vector_wrapper_tpetra_block<T>::value
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

//   auto Ctpb = *C.data(); //mv
//   auto Atpb = *(A.pressioObj()->data()); //v
//   auto Btpb = *B.data(); //mv

//   auto Ctp = Ctpb.getMultiVectorView();
//   using Atpb_t = mpl::remove_cvref_t<decltype(Atpb)>;
//   auto Atp = const_cast<Atpb_t &>(Atpb).getVectorView();
//   auto Btp = Btpb.getMultiVectorView();

//   Ctp.elementWiseMultiply(alpha, Atp, Btp, beta);
// }
#endif  // OPS_TPETRA_BLOCK_OPS_LEVEL3_HPP_
