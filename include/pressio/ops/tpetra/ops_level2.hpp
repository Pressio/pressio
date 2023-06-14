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

template <class A_type, class x_type, class y_type, class alpha_t, class beta_t>
::pressio::mpl::enable_if_t<
  ::pressio::is_multi_vector_tpetra<A_type>::value and
  ::pressio::is_vector_eigen<x_type>::value and
  ::pressio::is_vector_tpetra<y_type>::value
  >
_product_tpetra_mv_sharedmem_vec(const alpha_t & alpha,
				 const A_type & A,
				 const x_type & x,
				 const beta_t & beta,
				 y_type & y)
{
  assert(size_t(A.getNumVectors()) == size_t(::pressio::ops::extent(x,0)));

  using x_sc_t = typename ::pressio::Traits<x_type>::scalar_type;
  using kokkos_view_t = Kokkos::View<const x_sc_t*, Kokkos::HostSpace,
				     Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
  kokkos_view_t xview(x.data(), ::pressio::ops::extent(x,0));

  const auto ALocalView_h = A.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());
  const auto yLocalView_h = y.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  const char ctA = 'N';
  // Tpetra::Vector is implemented as a special case of MultiVector //
  // so getLocalView returns a rank-2 view so in order to get
  // view with rank==1 I need to explicitly get the subview of that
  const auto yLocalView_drank1 = Kokkos::subview(yLocalView_h, Kokkos::ALL(), 0);
  ::KokkosBlas::gemv(&ctA, alpha, ALocalView_h, xview, beta, yLocalView_drank1);
}

// when the operand is a kokkos wrapper we use kokkos functionalities directly
template <class A_type, class x_type, class y_type, class alpha_t, class beta_t>
::pressio::mpl::enable_if_t<
  ::pressio::is_multi_vector_tpetra<A_type>::value and
  ::pressio::is_vector_kokkos<x_type>::value and
  ::pressio::is_vector_tpetra<y_type>::value
  >
_product_tpetra_mv_sharedmem_vec(const alpha_t & alpha,
					const A_type & A,
					const x_type & x,
					const beta_t & beta,
					y_type & y)
{
  assert(x.span_is_contiguous());
  assert(size_t(A.getNumVectors()) == size_t(::pressio::ops::extent(x,0)));
  const char ctA = 'N';
  const auto ALocalView_d = A.getLocalViewDevice(Tpetra::Access::ReadOnlyStruct());

  // Tpetra::Vector is implemented as a special case of MultiVector //
  // so getLocalView returns a rank-2 view so in order to get
  // view with rank==1 I need to explicitly get the subview of that
  const auto yLocalView_drank2 = y.getLocalViewDevice(Tpetra::Access::ReadWriteStruct());
  const auto yLocalView_drank1 = Kokkos::subview(yLocalView_drank2, Kokkos::ALL(), 0);
  ::KokkosBlas::gemv(&ctA, alpha, ALocalView_d, x, beta, yLocalView_drank1);
}

template <class A_type, class x_type, class y_type, class alpha_t, class beta_t>
::pressio::mpl::enable_if_t<
  ::pressio::is_multi_vector_tpetra<A_type>::value and
  (::pressio::is_expression_acting_on_kokkos<x_type>::value or
   ::pressio::is_dense_vector_teuchos<x_type>::value) and
  ::pressio::is_vector_tpetra<y_type>::value
  >
_product_tpetra_mv_sharedmem_vec(const alpha_t & alpha,
					const A_type & A,
					const x_type & x,
					const beta_t & beta,
					y_type & y)
{
  assert(size_t(A.getNumVectors()) == size_t(::pressio::ops::extent(x, 0)));
  // const char ctA = 'N';
  const auto ALocalView_d = A.getLocalViewDevice(Tpetra::Access::ReadOnly);

  // Tpetra::Vector is implemented as a special case of MultiVector //
  // so getLocalView returns a rank-2 view so in order to get
  // view with rank==1 I need to explicitly get the subview of that
  const auto yLocalView_drank2 = y.getLocalViewDevice(Tpetra::Access::ReadWrite);
  const auto yLocalView_drank1 = Kokkos::subview(yLocalView_drank2, Kokkos::ALL(), 0);
  const auto x_size = ::pressio::ops::extent(x, 0);
  const auto y_size = ::pressio::ops::extent(yLocalView_drank1, 0);
  const auto zero = ::pressio::utils::Constants<typename pressio::Traits<x_type>::scalar_type>::zero();
  for (size_t i = 0; i < y_size; ++i) {
    if (beta == zero) {
      yLocalView_drank1(i) = zero;
    }
    else {
      yLocalView_drank1(i) *= beta;
    }
    if (alpha != zero) {
      for (size_t j = 0; j < x_size; ++j) {
        yLocalView_drank1(i) += alpha * ALocalView_d(i, j) * x(j);
      }
    }
  }
}
}//end namespace pressio::ops::impl


// -------------------------------
// y = beta * y + alpha*A*x
//
// x is Teuchos Vector
// A = tpetra::MultiVector
// y = tpetra vector
// -------------------------------
template < class A_type, class x_type, class y_type, class alpha_t, class beta_t>
::pressio::mpl::enable_if_t<
  // level2 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<x_type>::rank == 1
  && ::pressio::Traits<y_type>::rank == 1
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra<A_type>::value
  && ::pressio::is_vector_tpetra<y_type>::value
  && ::pressio::is_dense_vector_teuchos<x_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, x_type, y_type>::value
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  >
product(::pressio::nontranspose /*unused*/,
	const alpha_t & alpha,
	const A_type & A,
	const x_type & x,
	const beta_t & beta,
	y_type & y)
{
  ::pressio::ops::impl::_product_tpetra_mv_sharedmem_vec(alpha, A, x, beta, y);
}

// -------------------------------
// y = beta * y + alpha*A*x
//
// x is Kokkos Vector
// A = tpetra::MultiVector
// y = tpetra vector
// -------------------------------
template < class A_type, class x_type, class y_type, class alpha_t, class beta_t>
::pressio::mpl::enable_if_t<
  // level2 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<x_type>::rank == 1
  && ::pressio::Traits<y_type>::rank == 1
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra<A_type>::value
  && ::pressio::is_vector_tpetra<y_type>::value
  && (::pressio::is_vector_kokkos<x_type>::value
   || ::pressio::is_expression_acting_on_kokkos<x_type>::value)
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, x_type, y_type>::value
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  >
product(::pressio::nontranspose /*unused*/,
	const alpha_t & alpha,
	const A_type & A,
	const x_type & x,
	const beta_t & beta,
	y_type & y)
{
  ::pressio::ops::impl::_product_tpetra_mv_sharedmem_vec(alpha, A, x, beta, y);
}

// -------------------------------
// y = beta * y + alpha*A*x
//
// x is Eigen Vector
// A = tpetra::MultiVector
// y = tpetra vector
// -------------------------------
#ifdef PRESSIO_ENABLE_TPL_EIGEN
template < class A_type, class x_type, class y_type, class alpha_t, class beta_t>
::pressio::mpl::enable_if_t<
  // level2 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<x_type>::rank == 1
  && ::pressio::Traits<y_type>::rank == 1
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra<A_type>::value
  && ::pressio::is_vector_tpetra<y_type>::value
  && ::pressio::is_vector_eigen<x_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, x_type, y_type>::value
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  >
product(::pressio::nontranspose /*unused*/,
	const alpha_t & alpha,
	const A_type & A,
	const x_type & x,
	const beta_t & beta,
	y_type & y)
{

  //makesure x is contiguous
  assert(x.innerSize() == x.outerStride());
  ::pressio::ops::impl::_product_tpetra_mv_sharedmem_vec(alpha, A, x, beta, y);
}

// -------------------------------
// y = alpha*A*x, construct y
//
// x is Eigen Vector
// A = tpetra::MultiVector
// y = tpetra vector
// -------------------------------
template <class y_type, class A_type, class x_type, class alpha_t>
::pressio::mpl::enable_if_t<
  // level2 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<x_type>::rank == 1
  && ::pressio::Traits<y_type>::rank == 1
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra<A_type>::value
  && ::pressio::is_vector_tpetra<y_type>::value
  && ::pressio::is_vector_eigen<x_type>::value
// scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, x_type, y_type>::value
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value),
  y_type
  >
product(::pressio::nontranspose mode,
	const alpha_t & alpha,
	const A_type & A,
	const x_type & x)
{

  auto rowMap = A.getMap();
  y_type y(rowMap);
  using y_sc_t = typename y_type::scalar_type;
  product(mode, alpha, A, x, y_sc_t(0), y);
  return y;
}

// -------------------------------
// y = beta * y + alpha*A*x
//
// x is Pressio expression based on Eigen
// A = tpetra::MultiVector
// y = tpetra vector
// -------------------------------
template < class A_type, class x_type, class y_type, class alpha_t, class beta_t>
::pressio::mpl::enable_if_t<
  // level2 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<x_type>::rank == 1
  && ::pressio::Traits<y_type>::rank == 1
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra<A_type>::value
  && ::pressio::is_vector_tpetra<y_type>::value
  && ::pressio::is_expression_acting_on_eigen<x_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, x_type, y_type>::value
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  >
product(::pressio::nontranspose /*unused*/,
	const alpha_t & alpha,
	const A_type & A,
	const x_type & x,
	const beta_t & beta,
	y_type & y)
{
  const auto A_h = A.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());
  const auto y2_h = y.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
  const auto y_h = Kokkos::subview(y2_h, Kokkos::ALL(), 0);

  assert( ::pressio::ops::extent(y_h, 0) == ::pressio::ops::extent(A_h, 0) );
  assert( ::pressio::ops::extent(x, 0) == ::pressio::ops::extent(A_h, 1) );

  const auto zero = ::pressio::utils::Constants<beta_t>::zero();
  ::pressio::ops::scale(y_h, beta);
  if (alpha == zero) {
    return;
  }

  using sc_t = typename ::pressio::Traits<A_type>::scalar_type;
  const std::size_t m = ::pressio::ops::extent(A_h, 0);
  const std::size_t n = ::pressio::ops::extent(A_h, 1);
  Kokkos::parallel_for(m, KOKKOS_LAMBDA (const auto &i) {
    sc_t t{};
    for (std::size_t j = 0; j < n; j++) {
      t += alpha * A_h(i, j) * x(j);
    }
    y_h(i) = t;
  });
}

// -------------------------------
// y = beta * y + alpha*A^T*x
//
// x = tpetra Vector
// A = tpetra::MultiVector
// y = Eigen vector
// -------------------------------
template <class A_type, class x_type, class y_type, class alpha_t, class beta_t>
::pressio::mpl::enable_if_t<
  // level2 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<x_type>::rank == 1
  && ::pressio::Traits<y_type>::rank == 1
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra<A_type>::value
  && ::pressio::is_vector_tpetra<x_type>::value
  && (::pressio::is_vector_eigen<y_type>::value
   || ::pressio::is_expression_acting_on_eigen<y_type>::value)
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, x_type, y_type>::value
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  >
product(::pressio::transpose /*unused*/,
	const alpha_t & alpha,
	const A_type & A,
	const x_type & x,
	const beta_t & beta,
	y_type & y)
{
  // // dot product of each vector in A with vecB
  //  Apparently, trilinos does not support this...
  //    the native dot() method of multivectors is only for
  //    dot product of two multivectors with same # of columns.
  //    So we have to extract each column vector
  //    from A and do dot product one a time

  assert(size_t(A.getNumVectors()) == size_t(::pressio::ops::extent(y,0)));

  const auto zero = ::pressio::utils::Constants<beta_t>::zero();
  const auto numVecs = ::pressio::ops::extent(A, 1);
  for (std::size_t i=0; i<(std::size_t)numVecs; i++)
    {
      // colI is a Teuchos::RCP<Vector<...>>
      const auto colI = A.getVector(i);
      y(i) = beta == zero ? zero : beta * y(i);
      if (!(alpha == zero)) {
        y(i) += alpha * colI->dot(x);
      }
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
template <class A_type, class x_type, class y_type, class alpha_t, class beta_t>
::pressio::mpl::enable_if_t<
  // level2 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<x_type>::rank == 1
  && ::pressio::Traits<y_type>::rank == 1
  // TPL/container specific
  && ::pressio::is_multi_vector_tpetra<A_type>::value
  && ::pressio::is_vector_tpetra<x_type>::value
  && (::pressio::is_vector_kokkos<y_type>::value
   || ::pressio::is_expression_acting_on_kokkos<y_type>::value)
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, x_type, y_type>::value
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  >
product(::pressio::transpose /*unused*/,
	const alpha_t & alpha,
	const A_type & A,
	const x_type & x,
	const beta_t & beta,
	y_type & y)
{
  auto y_native = ::pressio::ops::impl::get_native(y);
  using y_native_type = typename ::pressio::mpl::remove_cvref<decltype(y_native)>::type;
  Kokkos::View<
    typename ::pressio::Traits<y_type>::scalar_type*,
    typename y_native_type::device_type
    > ATx("ATx", y.extent(0));
  auto request = Tpetra::idot(ATx, A, x);
  request->wait();
  ::KokkosBlas::axpby(alpha, ATx, beta, y_native);
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
