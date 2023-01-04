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

#ifndef OPS_TPETRA_BLOCK_OPS_LEVEL2_HPP_
#define OPS_TPETRA_BLOCK_OPS_LEVEL2_HPP_

#include "Tpetra_idot.hpp"
#include <KokkosBlas1_axpby.hpp>
#include "KokkosBlas2_gemv.hpp"

namespace pressio{ namespace ops{

// -------------------------------
// y = beta * y + alpha*A*x
//
// x is Teuchos Vector
// A = tpetra block MultiVector
// y = tpetra block vector
// -------------------------------
template < class A_type, class x_type, class y_type, class alpha_t, class beta_t>
::pressio::mpl::enable_if_t<
  // level2 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<x_type>::rank == 1
  && ::pressio::Traits<y_type>::rank == 1
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra_block<A_type>::value
  && ::pressio::is_vector_tpetra_block<y_type>::value
  && ::pressio::is_dense_vector_teuchos<x_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, x_type, y_type>::value
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  >
product(::pressio::nontranspose mode,
	const alpha_t & alpha,
	const A_type & A,
	const x_type & x,
	const beta_t & beta,
	y_type & y)
{

  auto A_tpetra_mv = A.getMultiVectorView();
  auto y_tpetra_v  = y.getVectorView();
  product(mode, alpha, A_tpetra_mv, x, beta, y_tpetra_v);
}

// -------------------------------
// y = beta * y + alpha*A*x
//
// x is Kokkos Vector
// A = tpetra block MultiVector
// y = tpetra block vector
// -------------------------------
template < class A_type, class x_type, class y_type, class alpha_t, class beta_t>
::pressio::mpl::enable_if_t<
  // level2 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<x_type>::rank == 1
  && ::pressio::Traits<y_type>::rank == 1
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra_block<A_type>::value
  && ::pressio::is_vector_tpetra_block<y_type>::value
  && (::pressio::is_vector_kokkos<x_type>::value
   || ::pressio::is_expression_acting_on_kokkos<x_type>::value)
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, x_type, y_type>::value
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  >
product(::pressio::nontranspose mode,
	const alpha_t & alpha,
	const A_type & A,
	const x_type & x,
	const beta_t & beta,
	y_type & y)
{

  auto A_tpetra_mv = A.getMultiVectorView();
  auto y_tpetra_v  = y.getVectorView();
  product(mode, alpha, A_tpetra_mv, x, beta, y_tpetra_v);
}

// -------------------------------
// y = beta * y + alpha*A*x
//
// x is Eigen Vector or Pressio expression based on Eigen container
// A = tpetra block MultiVector
// y = tpetra block vector
// -------------------------------
#ifdef PRESSIO_ENABLE_TPL_EIGEN
template < class A_type, class x_type, class y_type, class alpha_t, class beta_t>
::pressio::mpl::enable_if_t<
  // level2 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<x_type>::rank == 1
  && ::pressio::Traits<y_type>::rank == 1
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra_block<A_type>::value
  && ::pressio::is_vector_tpetra_block<y_type>::value
  && (::pressio::is_vector_eigen<x_type>::value
   || ::pressio::is_expression_acting_on_eigen<x_type>::value)
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, x_type, y_type>::value
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  >
product(::pressio::nontranspose mode,
	const alpha_t & alpha,
	const A_type & A,
	const x_type & x,
	const beta_t & beta,
	y_type & y)
{

  auto A_tpetra_mv = A.getMultiVectorView();
  auto y_tpetra_v  = y.getVectorView();
  product(mode, alpha, A_tpetra_mv, x, beta, y_tpetra_v);
}

// -------------------------------
// y = beta * y + alpha*A^T*x
//
// x = tpetra block Vector
// A = tpetra block MultiVector
// y = Eigen Vector or Pressio expression based on Eigen container
// -------------------------------
template < class A_type, class x_type, class y_type, class alpha_t, class beta_t>
::pressio::mpl::enable_if_t<
  // level2 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<x_type>::rank == 1
  && ::pressio::Traits<y_type>::rank == 1
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra_block<A_type>::value
  && ::pressio::is_vector_tpetra_block<x_type>::value
  && (::pressio::is_vector_eigen<y_type>::value
   || ::pressio::is_expression_acting_on_eigen<y_type>::value)
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, x_type, y_type>::value
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  >
product(::pressio::transpose mode,
	const alpha_t & alpha,
	const A_type & A,
	const x_type & x,
	const beta_t & beta,
	y_type & y)
{

  auto A_tpetra_mv = A.getMultiVectorView();
  auto x_tpetra_v  = static_cast<x_type>(x).getVectorView();
  product(mode, alpha, A_tpetra_mv, x_tpetra_v, beta, y);
}
#endif

// -------------------------------
// y = beta * y + alpha*A^T*x
//
// x = tpetra block Vector
// A = tpetra block MultiVector
// y = Kokkos vector
// -------------------------------
template < class A_type, class x_type, class y_type, class alpha_t, class beta_t>
::pressio::mpl::enable_if_t<
  // level2 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<x_type>::rank == 1
  && ::pressio::Traits<y_type>::rank == 1
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra_block<A_type>::value
  && ::pressio::is_vector_tpetra_block<x_type>::value
  && (::pressio::is_vector_kokkos<y_type>::value
   || ::pressio::is_expression_acting_on_kokkos<y_type>::value)
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, x_type, y_type>::value
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  >
product(::pressio::transpose mode,
	const alpha_t & alpha,
	const A_type & A,
	const x_type & x,
	const beta_t & beta,
	y_type & y)
{

  auto A_tpetra_mv = A.getMultiVectorView();
  auto x_tpetra_v  = static_cast<x_type>(x).getVectorView();
  product(mode, alpha, A_tpetra_mv, x_tpetra_v, beta, y);
}

}}//end namespace pressio::ops






// /* -------------------------------------------------------------------
//  * op(A) = A
//  * x is a sharedmem but NOT kokkos
//  *-------------------------------------------------------------------*/
// template < typename A_type, typename x_type, typename scalar_type, typename y_type>
// ::pressio::mpl::enable_if_t<
//   ::pressio::is_multi_vector_wrapper_tpetra_block<A_type>::value
//   and
//   (x_type::traits::wrapped_package_identifier != ::pressio::containers::details::WrappedPackageIdentifier::Kokkos) and
//   (::pressio::containers::predicates::is_vector_wrapper_eigen<x_type>::value or
//    ::pressio::containers::predicates::is_vector_wrapper_teuchos<x_type>::value or
//    ::pressio::containers::predicates::span_expression<x_type>::value)
//   and
//   ::pressio::containers::predicates::is_vector_wrapper_tpetra_block<y_type>::value
//   >
// product(::pressio::nontranspose mode,
// 	const scalar_type alpha,
// 	const A_type & A,
// 	const x_type & x,
// 	const scalar_type beta,
// 	y_type & y)
// {
//   using kokkos_view_t = Kokkos::View<const scalar_type*, Kokkos::HostSpace,
// 				     Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
//   kokkos_view_t xview(x.data()->data(), x.extent(0));

//   const auto ALocalView_h = A.data()->getMultiVectorView().getLocalViewHost();
//   const auto yLocalView_h = y.data()->getVectorView().getLocalViewHost();
//   const char ctA = 'N';
//   // Tpetra::Vector is implemented as a special case of MultiVector //
//   // so getLocalView returns a rank-2 view so in order to get
//   // view with rank==1 I need to explicitly get the subview of that
//   const auto yLocalView_drank1 = Kokkos::subview(yLocalView_h, Kokkos::ALL(), 0);
//   ::KokkosBlas::gemv(&ctA, alpha, ALocalView_h, xview, beta, yLocalView_drank1);
// }

// /* -------------------------------------------------------------------
//  * x is a distributed Tpetra block vector wrapper
//  *-------------------------------------------------------------------*/
// // y = sharedmem vec not kokkos
// template <typename A_type, typename x_type, typename y_type, typename scalar_type>
// ::pressio::mpl::enable_if_t<
//   ::pressio::containers::predicates::is_multi_vector_wrapper_tpetra_block<A_type>::value
//   and
//   ::pressio::containers::predicates::is_vector_wrapper_tpetra_block<x_type>::value
//   and
//   ::pressio::ops::constraints::sharedmem_host_subscriptable_rank1_container<y_type>::value and
//   (y_type::traits::wrapped_package_identifier !=
//    ::pressio::containers::details::WrappedPackageIdentifier::Kokkos)
//   >
// product(::pressio::transpose mode,
// 	const scalar_type alpha,
// 	const A_type & A,
// 	const x_type & x,
// 	const scalar_type beta,
// 	y_type & y)
// {

//   /* workaround the non-constness of getVectorView*/
//   using wrapped_t = typename ::pressio::containers::details::traits<x_type>::wrapped_t;
//   using ord_t = typename ::pressio::containers::details::traits<A_type>::global_ordinal_t;
//   const auto xvv = const_cast<wrapped_t*>(x.data())->getVectorView();
//   const auto mvA_mvv = A.data()->getMultiVectorView();
//   const auto numVecs = A.numVectors();
//   for (ord_t i=0; i<numVecs; i++){
//     // colI is a Teuchos::RCP<Vector<...>>
//     const auto colI = mvA_mvv.getVector(i);
//     y(i) = beta*y(i) + alpha * colI->dot(xvv);
//   }
// }

// // y = wrapper of Kokkos vector
// template <typename A_type, typename x_type, typename y_type, typename scalar_type>
// ::pressio::mpl::enable_if_t<
//   ::pressio::containers::predicates::is_multi_vector_wrapper_tpetra_block<A_type>::value
//   and
//   ::pressio::containers::predicates::is_vector_wrapper_tpetra_block<x_type>::value
//   and
//   ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<y_type>::value
//   >
// product(::pressio::transpose mode,
// 	const scalar_type alpha,
// 	const A_type & A,
// 	const x_type & x,
// 	const scalar_type beta,
// 	y_type & y)
// {
//   static_assert
//     (std::is_same<
//      typename ::pressio::containers::details::traits<A_type>::device_t,
//      typename ::pressio::containers::details::traits<x_type>::device_t>::value,
//      "Tpetra MV dot V: operands do not have the same device type");

//   static_assert
//     (std::is_same<
//      typename ::pressio::containers::details::traits<x_type>::device_t,
//      typename ::pressio::containers::details::traits<y_type>::device_t>::value,
//      "Tpetra MV dot V: V and result do not have the same device type");

//   using kokkos_v_t = typename ::pressio::containers::details::traits<y_type>::wrapped_t;
//   using v_t = ::pressio::containers::Vector<kokkos_v_t>;
//   using tpetra_blockvector_t = typename ::pressio::containers::details::traits<x_type>::wrapped_t;

//   const auto A_mvv = A.data()->getMultiVectorView();
//   const auto x_vv = const_cast<tpetra_blockvector_t*>(x.data())->getVectorView();

//   v_t ATx(y.extent(0));
//   auto request = Tpetra::idot(*ATx.data(), A_mvv, x_vv);
//   request->wait();
//   ::KokkosBlas::axpby(alpha, *ATx.data(), beta, *y.data());
// }
#endif  // OPS_TPETRA_BLOCK_OPS_LEVEL2_HPP_
