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
template <typename A_type, typename B_type, typename ScalarType, typename C_type>
::pressio::mpl::enable_if_t<
  ::pressio::is_dense_matrix_eigen<A_type>::value and
  ::pressio::is_dense_matrix_eigen<B_type>::value and
  ::pressio::is_dense_matrix_eigen<C_type>::value
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const ScalarType alpha,
	const A_type & A,
	const B_type & B,
	const ScalarType beta,
	C_type & C)
{
  static_assert
    (::pressio::are_scalar_compatible<A_type, B_type, C_type>::value,
     "Types are not scalar compatible");

  assert( ::pressio::ops::extent(C, 0) == ::pressio::ops::extent(A, 1) );
  assert( ::pressio::ops::extent(C, 1) == ::pressio::ops::extent(B, 1) );
  assert( ::pressio::ops::extent(A, 0) == ::pressio::ops::extent(B, 0) );
  C = beta * C + alpha * A.transpose() * B;
}

// template < typename A_type, typename B_type, typename ScalarType, typename C_type>
// ::pressio::mpl::enable_if_t<
//   ::pressio::is_dense_matrix_eigen<B_type>::value and
//   ::pressio::is_dense_matrix_eigen<C_type>::value
//   >
// product(::pressio::transpose,
// 	::pressio::nontranspose,
// 	const ScalarType alpha,
// 	const ::pressio::containers::experimental::MultiVectorSet<A_type> & A,
// 	const B_type & B,
// 	const ScalarType beta,
// 	C_type & C)
// {
//   static_assert
//     (::pressio::are_scalar_compatible
//      <::pressio::containers::experimental::MultiVectorSet<A_type>, B_type, C_type>::value,
//      "Types are not scalar compatible");

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

//-------------------------------------------
// specialize for op(A) = A and op(B) = B
//-------------------------------------------
template <typename A_type, typename B_type, typename ScalarType, typename C_type>
::pressio::mpl::enable_if_t<
  ::pressio::is_dense_matrix_eigen<A_type>::value and
  ::pressio::is_dense_matrix_eigen<B_type>::value and
  ::pressio::is_dense_matrix_eigen<C_type>::value
  >
product(::pressio::nontranspose modeA,
	::pressio::nontranspose modeB,
	const ScalarType alpha,
	const A_type & A,
	const B_type & B,
	const ScalarType beta,
	C_type & C)
{
  static_assert
    (::pressio::are_scalar_compatible<A_type, B_type, C_type>::value,
     "Types are not scalar compatible");

  assert( ::pressio::ops::extent(C, 0) == ::pressio::ops::extent(A, 0) );
  assert( ::pressio::ops::extent(C, 1) == ::pressio::ops::extent(B, 1) );
  assert( ::pressio::ops::extent(A, 1) == ::pressio::ops::extent(B, 0) );
  C = beta * C + alpha * A * B;
}

// template < typename A_type, typename B_type, typename ScalarType, typename C_type>
// ::pressio::mpl::enable_if_t<
//   ::pressio::is_dense_matrix_eigen<B_type>::value and
//   ::pressio::is_dense_matrix_eigen<C_type>::value
//   >
// product(::pressio::nontranspose,
// 	::pressio::nontranspose,
// 	const ScalarType alpha,
// 	const ::pressio::containers::experimental::MultiVectorSet<A_type> & A,
// 	const B_type & B,
// 	const ScalarType beta,
// 	C_type & C)
// {
//   static_assert
//     (::pressio::are_scalar_compatible
//      <::pressio::containers::experimental::MultiVectorSet<A_type>, B_type, C_type>::value,
//      "Types are not scalar compatible");

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

/***********************************
* special case A==B and op(A) = transpose
**********************************/
template <typename A_type, typename ScalarType, typename C_type>
::pressio::mpl::enable_if_t<
  ::pressio::is_dense_matrix_eigen<A_type>::value and
  ::pressio::is_dense_matrix_eigen<C_type>::value
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const ScalarType alpha,
	const A_type & A,
	const ScalarType beta,
	C_type & C)
{
  static_assert
    (::pressio::are_scalar_compatible<A_type, C_type>::value,
     "Types are not scalar compatible");

  // const auto & AE = *A.data();
  // auto & CE = *C.data();
  C = beta * C + alpha * A.transpose() * A;
}

template <typename C_type, typename A_type, typename ScalarType>
::pressio::mpl::enable_if_t<
  ::pressio::is_dense_matrix_eigen<A_type>::value and
  ::pressio::is_dense_matrix_eigen<C_type>::value,
  C_type
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const ScalarType alpha,
	const A_type & A)
{
  static_assert
    (::pressio::are_scalar_compatible<A_type, C_type>::value,
     "Types are not scalar compatible");

  constexpr auto zero = ::pressio::utils::constants<ScalarType>::zero();
  C_type C(::pressio::ops::extent(A, 1), ::pressio::ops::extent(A, 1));
  product(modeA, modeB, alpha, A, A, zero, C);
  return C;
}

//-------------------------------------------
// C = beta * C + alpha*A*B
// specialize for when A = asDiagonalMatrix expression
//-------------------------------------------
template <typename A_type, typename B_type, typename ScalarType, typename C_type>
::pressio::mpl::enable_if_t<
  ::pressio::is_dense_matrix_eigen<B_type>::value and
  ::pressio::is_dense_matrix_eigen<C_type>::value and 
  ::pressio::is_expression_asdiagonalmatrix<A_type>::value and 
  ::pressio::traits<A_type>::package_identifier == PackageIdentifier::Eigen
  >
product(::pressio::nontranspose modeA,
	::pressio::nontranspose modeB,
	const ScalarType alpha,
	const A_type & A,
	const B_type & B,
	const ScalarType beta,
	C_type & C)
{
  static_assert
    (::pressio::are_scalar_compatible<A_type, B_type, C_type>::value,
     "Types are not scalar compatible");

  assert( ::pressio::ops::extent(C, 0) == ::pressio::ops::extent(A, 0) );
  assert( ::pressio::ops::extent(C, 1) == ::pressio::ops::extent(B, 1) );
  assert( ::pressio::ops::extent(A, 1) == ::pressio::ops::extent(B, 0) );
  const auto & AE = *A.pressioObj();
  C = beta*C + alpha * (AE.asDiagonal() * B);
}

}}//end namespace pressio::ops
#endif  // OPS_EIGEN_OPS_LEVEL3_HPP_
