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

#ifndef OPS_TPETRA_OPS_LEVEL2_HPP_
#define OPS_TPETRA_OPS_LEVEL2_HPP_

#include "Tpetra_idot.hpp"
#include <KokkosBlas1_axpby.hpp>
#include "KokkosBlas2_gemv.hpp"

namespace pressio{ namespace ops{


/*
y = beta * y + alpha*op(A)*x

op(A) = A or A^T
*/

namespace impl{

template <typename A_type, typename x_type, typename y_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  ::pressio::is_multi_vector_tpetra<A_type>::value and
  ::pressio::is_vector_tpetra<y_type>::value
  >
_product_tpetra_mv_sharedmem_vec(const scalar_type alpha,
				 const A_type & A,
				 const x_type & x,
				 const scalar_type beta,
				 y_type & y)
{
  assert(size_t(A.getNumVectors()) == size_t(::pressio::ops::extent(x,0)));

  using kokkos_view_t = Kokkos::View<const scalar_type*, Kokkos::HostSpace,
				     Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
  kokkos_view_t xview(x.data(), ::pressio::ops::extent(x,0));

  const auto ALocalView_h = A.getLocalViewHost();
  const auto yLocalView_h = y.getLocalViewHost();
  const char ctA = 'N';
  // Tpetra::Vector is implemented as a special case of MultiVector //
  // so getLocalView returns a rank-2 view so in order to get
  // view with rank==1 I need to explicitly get the subview of that
  const auto yLocalView_drank1 = Kokkos::subview(yLocalView_h, Kokkos::ALL(), 0);
  ::KokkosBlas::gemv(&ctA, alpha, ALocalView_h, xview, beta, yLocalView_drank1);
}

// when the operand is a kokkos wrapper we use kokkos functionalities directly
template <typename A_type, typename x_type, typename y_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  ::pressio::is_multi_vector_tpetra<A_type>::value and
  ::pressio::is_vector_tpetra<y_type>::value
  >
_product_tpetra_mv_sharedmem_vec_kokkos(const scalar_type alpha,
					const A_type & A,
					const x_type & x,
					const scalar_type beta,
					y_type & y)
{
  // make sure the tpetra mv has same exe space of the kokkos vector wrapper
  using tpetra_mv_dev_t = typename ::pressio::Traits<A_type>::device_type;
  using kokkos_v_dev_t  = typename ::pressio::Traits<x_type>::device_type;
  static_assert
    ( std::is_same<tpetra_mv_dev_t, kokkos_v_dev_t>::value,
      "product: tpetra MV and kokkos wrapper need to have same device type" );

  assert(size_t(A.getNumVectors()) == size_t(::pressio::ops::extent(x,0)));
  const char ctA = 'N';
  const auto ALocalView_d = A.getLocalViewDevice();

  // Tpetra::Vector is implemented as a special case of MultiVector //
  // so getLocalView returns a rank-2 view so in order to get
  // view with rank==1 I need to explicitly get the subview of that
  const auto yLocalView_drank2 = y.getLocalViewDevice();
  const auto yLocalView_drank1 = Kokkos::subview(yLocalView_drank2, Kokkos::ALL(), 0);
  ::KokkosBlas::gemv(&ctA, alpha, ALocalView_d, x, beta, yLocalView_drank1);
}
}//end namespace pressio::ops::impl


// -------------------------------
// y = beta * y + alpha*A*x
//
// x is Teuchos Vector
// A = tpetra::MultiVector
// y = tpetra vector
// -------------------------------
template < typename A_type, typename x_type, typename scalar_type, typename y_type>
::pressio::mpl::enable_if_t<
  ::pressio::is_multi_vector_tpetra<A_type>::value
  and ::pressio::is_vector_tpetra<y_type>::value
  and ::pressio::is_dense_vector_teuchos<x_type>::value
  >
product(::pressio::nontranspose mode,
	const scalar_type alpha,
	const A_type & A,
	const x_type & x,
	const scalar_type beta,
	y_type & y)
{
  static_assert
    (::pressio::are_scalar_compatible<A_type, x_type, y_type>::value,
     "Types are not scalar compatible");

  ::pressio::ops::impl::_product_tpetra_mv_sharedmem_vec(alpha, A, x, beta, y);
}

// -------------------------------
// y = beta * y + alpha*A*x
//
// x is Kokkos Vector
// A = tpetra::MultiVector
// y = tpetra vector
// -------------------------------
template < typename A_type, typename x_type, typename scalar_type, typename y_type>
::pressio::mpl::enable_if_t<
  ::pressio::is_multi_vector_tpetra<A_type>::value and
  ::pressio::is_vector_tpetra<y_type>::value and
  ::pressio::is_vector_kokkos<x_type>::value
  >
product(::pressio::nontranspose,
	const scalar_type alpha,
	const A_type & A,
	const x_type & x,
	const scalar_type beta,
	y_type & y)
{
  static_assert
    (::pressio::are_scalar_compatible<A_type, x_type, y_type>::value,
     "Types are not scalar compatible");
  assert(x.span_is_contiguous());

  ::pressio::ops::impl::_product_tpetra_mv_sharedmem_vec_kokkos(alpha, A, x, beta, y);
}

// -------------------------------
// y = beta * y + alpha*A*x
//
// x is Eigen Vector
// A = tpetra::MultiVector
// y = tpetra vector
// -------------------------------
#ifdef PRESSIO_ENABLE_TPL_EIGEN
template < typename A_type, typename x_type, typename scalar_type, typename y_type>
::pressio::mpl::enable_if_t<
  ::pressio::is_multi_vector_tpetra<A_type>::value
  and ::pressio::is_vector_tpetra<y_type>::value
  and ::pressio::is_vector_eigen<x_type>::value
  >
product(::pressio::nontranspose mode,
	const scalar_type alpha,
	const A_type & A,
	const x_type & x,
	const scalar_type beta,
	y_type & y)
{
  static_assert
    (::pressio::are_scalar_compatible<A_type, x_type, y_type>::value,
     "Types are not scalar compatible");

  //makesure x is contiguous
  assert(x.innerSize() == x.outerStride());
  ::pressio::ops::impl::_product_tpetra_mv_sharedmem_vec(alpha, A, x, beta, y);
}

// -------------------------------
// y = beta * y + alpha*A^T*x
//
// x = tpetra Vector
// A = tpetra::MultiVector
// y = Eigen vector
// -------------------------------
template <typename A_type, typename x_type, typename y_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  ::pressio::is_multi_vector_tpetra<A_type>::value
  and ::pressio::is_vector_tpetra<x_type>::value
  and ::pressio::is_vector_eigen<y_type>::value
  >
product(::pressio::transpose mode,
	const scalar_type alpha,
	const A_type & A,
	const x_type & x,
	const scalar_type beta,
	y_type & y)
{
  // // dot product of each vector in A with vecB
  //  Apparently, trilinos does not support this...
  //    the native dot() method of multivectors is only for
  //    dot product of two multivectors with same # of columns.
  //    So we have to extract each column vector
  //    from A and do dot product one a time

  static_assert
    (::pressio::are_scalar_compatible<A_type, x_type, y_type>::value,
     "Types are not scalar compatible");

  assert(size_t(A.getNumVectors()) == size_t(::pressio::ops::extent(y,0)));

  const auto numVecs = ::pressio::ops::extent(A, 1);
  for (std::size_t i=0; i<(std::size_t)numVecs; i++)
    {
      // colI is a Teuchos::RCP<Vector<...>>
      const auto colI = A.getVector(i);
      y(i) = beta * y(i) + alpha * colI->dot(x);
    }
}
#endif


// -------------------------------
// y = beta * y + alpha*A^T*x
//
// x = tpetra Vector
// A = tpetra::MultiVector
// y = Kokkos vector
// -------------------------------
template <typename A_type, typename x_type, typename y_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  ::pressio::is_multi_vector_tpetra<A_type>::value and
  ::pressio::is_vector_tpetra<x_type>::value and
  ::pressio::is_vector_kokkos<y_type>::value
  >
product(::pressio::transpose mode,
	const scalar_type alpha,
	const A_type & A,
	const x_type & x,
	const scalar_type beta,
	y_type & y)
{
  static_assert
    (::pressio::are_scalar_compatible<A_type, x_type, y_type>::value,
     "Tpetra MV dot V: operands do not have matching scalar type");

  static_assert
    (std::is_same<
     typename ::pressio::Traits<A_type>::device_type,
     typename ::pressio::Traits<x_type>::device_type>::value,
     "Tpetra MV dot V: operands do not have the same device type");

  static_assert
    (std::is_same<
     typename ::pressio::Traits<x_type>::device_type,
     typename ::pressio::Traits<y_type>::device_type>::value,
     "Tpetra MV dot V: V and result do not have the same device type");

  y_type ATx("ATx", y.extent(0));
  auto request = Tpetra::idot(ATx, A, x);
  request->wait();
  ::KokkosBlas::axpby(alpha, ATx, beta, y);
}

}}//end namespace pressio::ops





// #ifdef PRESSIO_ENABLE_TPL_EIGEN
// /* -------------------------------------------------------------------
//  * y = beta * y + alpha*A*x

//  * x is a sharedmem vector wrapper eigen
//  * A = Eigen::MultiVector wrapper
//  * y = Tpetra Vector
//  * this covers the case where the matrix A acts locally to x
//    while y is distributed, so A*x only fills the corresponding part of y
//  *-------------------------------------------------------------------*/
// template < typename A_type, typename x_type, typename y_type, typename scalar_type>
// ::pressio::mpl::enable_if_t<
//       ::pressio::is_multi_vector_eigen<A_type>::value
//   and ::pressio::is_vector_eigen<x_type>::value
//   and ::pressio::is_vector_tpetra<y_type>::value
//   >
// product(::pressio::nontranspose mode,
// 	const scalar_type alpha,
// 	const A_type & A,
// 	const x_type & x,
// 	const scalar_type beta,
// 	y_type & y)
// {
//   static_assert
//     (are_scalar_compatible<A_type, x_type, y_type>::value,
//      "Types are not scalar compatible");
//   static_assert
//     (mpl::is_same<scalar_type, typename ::pressio::Traits<x_type>::scalar_type>::value,
//      "Scalar compatibility broken");

//   using kuv = Kokkos::MemoryTraits<Kokkos::Unmanaged>;

//   using A_view_t = Kokkos::View<const scalar_type**, Kokkos::HostSpace, kuv>;
//   A_view_t Aview(A.data()->data(), A.extent(0), A.extent(1));

//   using x_view_t = Kokkos::View<const scalar_type*, Kokkos::HostSpace, kuv>;
//   x_view_t xview(x.data()->data(), x.extent(0));

//   const char ctA = 'N';
//   const auto yLocalView_h = y.data()->getLocalViewHost();
//   const auto yLocalView_rank1 = Kokkos::subview(yLocalView_h, Kokkos::ALL(), 0);
//   ::KokkosBlas::gemv(&ctA, alpha, Aview, xview, beta, yLocalView_rank1);
// }
// #endif

// #ifdef PRESSIO_ENABLE_TPL_EIGEN
// /* -------------------------------------------------------------------
//  * y = beta * y + alpha*A^T*x
//  * x is a distributed Tpetra vector wrapper
//  * A is an Eigen MV
//  * y is vector wrapper Eigen
//  *-------------------------------------------------------------------*/
// template <typename A_type, typename x_type, typename y_type, typename scalar_type>
// ::pressio::mpl::enable_if_t<
//   ::pressio::is_vector_wpetra<x_type>::value and
//   ::pressio::is_multi_vector_eigen<A_type>::value and
//   ::pressio::is_vector_eigen<y_type>::value
//   >
// product(::pressio::transpose mode,
// 	const scalar_type alpha,
// 	const A_type & A,
// 	const x_type & x,
// 	const scalar_type beta,
// 	y_type & y)
// {
//   static_assert
//     (containers::predicates::are_scalar_compatible<A_type, x_type, y_type>::value,
//      "Types are not scalar compatible");
//   static_assert
//     (mpl::is_same<
//      scalar_type, typename ::pressio::containers::details::traits<x_type>::scalar_t>::value,
//      "Scalar compatibility broken");

//   using kuv = Kokkos::MemoryTraits<Kokkos::Unmanaged>;

//   using A_view_t = Kokkos::View<const scalar_type**, Kokkos::HostSpace, kuv>;
//   A_view_t Aview(A.data()->data(), A.extent(0), A.extent(1));

//   const auto xLocalView_h     = x.data()->getLocalViewHost();
//   const auto xLocalView_rank1 = Kokkos::subview(xLocalView_h, Kokkos::ALL(), 0);

//   using y_view_t = Kokkos::View<scalar_type*, Kokkos::HostSpace, kuv>;
//   y_view_t yview(y.data()->data(), y.extent(0));

//   const char ctA = 'T';
//   ::KokkosBlas::gemv(&ctA, alpha, Aview, xLocalView_rank1, beta, yview);
// }
// #endif
#endif  // OPS_TPETRA_OPS_LEVEL2_HPP_
