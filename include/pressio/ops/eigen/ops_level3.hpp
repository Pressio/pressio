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

#ifndef OPS_EIGEN_OPS_LEVEL3_HPP_
#define OPS_EIGEN_OPS_LEVEL3_HPP_

namespace pressio{ namespace ops{

/*
 * C = beta * C + alpha*op(A)*op(B)
*/

//-------------------------------------------
// specialize for op(A) = A^T and op(B) = B
//-------------------------------------------
template <
  class A_type, class B_type, class C_type,
  class alpha_t, class beta_t
  >
::pressio::mpl::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<B_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_native_container_eigen<A_type>::value
  && ::pressio::is_native_container_eigen<B_type>::value
  && ::pressio::is_native_container_eigen<C_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, B_type, C_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  /* constrained via is_convertible because the impl is using
     native Eigen expression whcih only work if the scalars are
     convertible to object scalar types*/
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

  assert( ::pressio::ops::extent(C, 0) == ::pressio::ops::extent(A, 1) );
  assert( ::pressio::ops::extent(C, 1) == ::pressio::ops::extent(B, 1) );
  assert( ::pressio::ops::extent(A, 0) == ::pressio::ops::extent(B, 0) );
  if (beta == pressio::utils::Constants<beta_t>::zero()) {
    C = alpha * A.transpose() * B;
  }
  else {
    C = beta * C + alpha * A.transpose() * B;
  }
}

//-------------------------------------------
// specialize for op(A) = A and op(B) = B
//-------------------------------------------
template <
  class A_type, class B_type, class C_type,
  class alpha_t, class beta_t
  >
::pressio::mpl::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<B_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_native_container_eigen<A_type>::value
  && ::pressio::is_native_container_eigen<B_type>::value
  && ::pressio::is_native_container_eigen<C_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, B_type, C_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  /* constrained via is_convertible because the impl is using
     native Eigen expression whcih only work if the scalars are
     convertible to object scalar types*/
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t,  typename ::pressio::Traits<A_type>::scalar_type>::value
  >
product(::pressio::nontranspose /*unused*/,
	::pressio::nontranspose /*unused*/,
	const alpha_t & alpha,
	const A_type & A,
	const B_type & B,
	const beta_t & beta,
	C_type & C)
{

  assert( ::pressio::ops::extent(C, 0) == ::pressio::ops::extent(A, 0) );
  assert( ::pressio::ops::extent(C, 1) == ::pressio::ops::extent(B, 1) );
  assert( ::pressio::ops::extent(A, 1) == ::pressio::ops::extent(B, 0) );
  if (beta == pressio::utils::Constants<beta_t>::zero()) {
    C = alpha * A * B;
  }
  else {
    C = beta * C + alpha * A * B;
  }
}

/***********************************
* special case A==B and op(A) = transpose
**********************************/
template <class A_type, class C_type, class alpha_t, class beta_t>
::pressio::mpl::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_native_container_eigen<A_type>::value
  && ::pressio::is_native_container_eigen<C_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, C_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  /* constrained via is_convertible because the impl is using
     native Eigen expression whcih only work if the scalars are
     convertible to object scalar types*/
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

  if (beta == pressio::utils::Constants<beta_t>::zero()) {
    C = alpha * A.transpose() * A;
  } else {
    C = beta * C + alpha * A.transpose() * A;
  }
}

template <class C_type, class A_type, class alpha_t>
::pressio::mpl::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_native_container_eigen<A_type>::value
  && ::pressio::is_native_container_eigen<C_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, C_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  /* constrained via is_convertible because the impl is using
     native Eigen expression whcih only work if the scalars are
     convertible to object scalar types*/
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value,
  C_type
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const alpha_t & alpha,
	const A_type & A)
{

  using sc_t = typename ::pressio::Traits<C_type>::scalar_type;
  constexpr auto zero = ::pressio::utils::Constants<sc_t>::zero();
  C_type C(::pressio::ops::extent(A, 1), ::pressio::ops::extent(A, 1));
  product(modeA, modeB, alpha, A, A, zero, C);
  return C;
}

//-------------------------------------------
// C = beta * C + alpha*A*B
// specialize for when A = asDiagonalMatrix expression
//-------------------------------------------
template <
  class A_type, class B_type, class C_type,
  class alpha_t, class beta_t
  >
::pressio::mpl::enable_if_t<
  // level3 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<B_type>::rank == 2
  && ::pressio::Traits<C_type>::rank == 2
  // TPL/container specific
  && ::pressio::is_expression_asdiagonalmatrix<A_type>::value
  && ::pressio::is_expression_acting_on_eigen<A_type>::value
  && ::pressio::is_native_container_eigen<B_type>::value
  && ::pressio::is_native_container_eigen<C_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, B_type, C_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  /* constrained via is_convertible because the impl is using
     native Eigen expression whcih only work if the scalars are
     convertible to object scalar types*/
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t,  typename ::pressio::Traits<A_type>::scalar_type>::value
  >
product(::pressio::nontranspose /*unused*/,
	::pressio::nontranspose /*unused*/,
	const alpha_t & alpha,
	const A_type & A,
	const B_type & B,
	const beta_t & beta,
	C_type & C)
{

  assert( ::pressio::ops::extent(C, 0) == ::pressio::ops::extent(A, 0) );
  assert( ::pressio::ops::extent(C, 1) == ::pressio::ops::extent(B, 1) );
  assert( ::pressio::ops::extent(A, 1) == ::pressio::ops::extent(B, 0) );
  if (beta == pressio::utils::Constants<beta_t>::zero()) {
    C = alpha * A.native() * B;
  } else {
    C = beta*C + alpha * A.native() * B;
  }
}

}}//end namespace pressio::ops




// template < class A_type, class B_type, class ScalarType, class C_type>
// ::pressio::mpl::enable_if_t<
//   ::pressio::is_dense_matrix_eigen<B_type>::value and
//   ::pressio::is_dense_matrix_eigen<C_type>::value
//   >
// product(::pressio::transpose,
//  ::pressio::nontranspose,
//  const ScalarType alpha,
//  const ::pressio::containers::experimental::MultiVectorSet<A_type> & A,
//  const B_type & B,
//  const ScalarType beta,
//  C_type & C)
// {

//   const auto & BE = *B.data();
//   auto & CE = *C.data();
//   for (std::size_t i=0; i<A.size(); ++i)
//   {
//     const auto currMatrixEigen = *(A(i).data());
//     assert( ::pressio::ops::extent(C, 0) == currMatrixEigen.cols() );
//     assert( ::pressio::ops::extent(B, 0) == currMatrixEigen.rows() );
//     CE.col(i) = beta * CE.col(i) + alpha * currMatrixEigen.transpose() * BE.col(i);
//   }
// }

// template < class A_type, class B_type, class ScalarType, class C_type>
// ::pressio::mpl::enable_if_t<
//   ::pressio::is_dense_matrix_eigen<B_type>::value and
//   ::pressio::is_dense_matrix_eigen<C_type>::value
//   >
// product(::pressio::nontranspose,
//  ::pressio::nontranspose,
//  const ScalarType alpha,
//  const ::pressio::containers::experimental::MultiVectorSet<A_type> & A,
//  const B_type & B,
//  const ScalarType beta,
//  C_type & C)
// {

//   const auto & BE = *B.data();
//   auto & CE = *C.data();
//   for (std::size_t i=0; i<A.size(); ++i)
//   {
//     const auto & currMatrixEigen = *(A(i).data());
//     assert( ::pressio::ops::extent(C, 0) == currMatrixEigen.rows() );
//     assert( ::pressio::ops::extent(B, 0) == currMatrixEigen.cols() );
//     CE.col(i) = beta * CE.col(i) + alpha * currMatrixEigen * BE.col(i);
//   }
// }
#endif  // OPS_EIGEN_OPS_LEVEL3_HPP_
