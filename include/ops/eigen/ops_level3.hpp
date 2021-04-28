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
 *
*/

//-------------------------------------------
// specialize for op(A) = A^T and op(B) = B
//-------------------------------------------
template <typename A_type, typename B_type, typename scalar_type, typename C_type>
::pressio::mpl::enable_if_t<
  ::pressio::ops::constraints::rank2_container_eigen_with_native_data_access<A_type>::value and
  ::pressio::ops::constraints::rank2_container_eigen_with_native_data_access<B_type>::value and
  ::pressio::ops::constraints::rank2_container_eigen_with_native_data_access<C_type>::value
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
    (containers::predicates::are_scalar_compatible<A_type, B_type, C_type>::value,
     "Types are not scalar compatible");

  assert( C.extent(0) == A.extent(1) );
  assert( C.extent(1) == B.extent(1) );
  assert( A.extent(0) == B.extent(0) );

  const auto & AE = *A.data();
  const auto & BE = *B.data();
  auto & CE = *C.data();
  CE = beta * CE + alpha * AE.transpose() * BE;
}

template < typename A_type, typename B_type, typename scalar_type, typename C_type>
::pressio::mpl::enable_if_t<
  ::pressio::ops::constraints::rank2_container_eigen_with_native_data_access<B_type>::value and
  ::pressio::ops::constraints::rank2_container_eigen_with_native_data_access<C_type>::value
  >
product(::pressio::transpose,
	::pressio::nontranspose,
	const scalar_type alpha,
	const ::pressio::containers::experimental::MultiVectorSet<A_type> & A,
	const B_type & B,
	const scalar_type beta,
	C_type & C)
{
  static_assert
    (containers::predicates::are_scalar_compatible
     <::pressio::containers::experimental::MultiVectorSet<A_type>, B_type, C_type>::value,
     "Types are not scalar compatible");

  const auto & BE = *B.data();
  auto & CE = *C.data();
  for (std::size_t i=0; i<A.size(); ++i)
  {
    const auto currMatrixEigen = *(A(i).data());
    assert( C.extent(0) == currMatrixEigen.cols() );
    assert( B.extent(0) == currMatrixEigen.rows() );
    CE.col(i) = beta * CE.col(i) + alpha * currMatrixEigen.transpose() * BE.col(i);
  }
}

//-------------------------------------------
// specialize for op(A) = A and op(B) = B
//-------------------------------------------
template <typename A_type, typename B_type, typename scalar_type, typename C_type>
::pressio::mpl::enable_if_t<
  ::pressio::ops::constraints::rank2_container_eigen_with_native_data_access<A_type>::value and
  ::pressio::ops::constraints::rank2_container_eigen_with_native_data_access<B_type>::value and
  ::pressio::ops::constraints::rank2_container_eigen_with_native_data_access<C_type>::value
  >
product(::pressio::nontranspose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const A_type & A,
	const B_type & B,
	const scalar_type beta,
	C_type & C)
{
  static_assert
    (containers::predicates::are_scalar_compatible<A_type, B_type, C_type>::value,
     "Types are not scalar compatible");

  assert( C.extent(0) == A.extent(0) );
  assert( C.extent(1) == B.extent(1) );
  assert( A.extent(1) == B.extent(0) );

  const auto & AE = *A.data();
  const auto & BE = *B.data();
  auto & CE = *C.data();

  CE = beta * CE + alpha * AE * BE;
}

template < typename A_type, typename B_type, typename scalar_type, typename C_type>
::pressio::mpl::enable_if_t<
  ::pressio::ops::constraints::rank2_container_eigen_with_native_data_access<B_type>::value and
  ::pressio::ops::constraints::rank2_container_eigen_with_native_data_access<C_type>::value
  >
product(::pressio::nontranspose,
	::pressio::nontranspose,
	const scalar_type alpha,
	const ::pressio::containers::experimental::MultiVectorSet<A_type> & A,
	const B_type & B,
	const scalar_type beta,
	C_type & C)
{
  static_assert
    (containers::predicates::are_scalar_compatible
     <::pressio::containers::experimental::MultiVectorSet<A_type>, B_type, C_type>::value,
     "Types are not scalar compatible");

  const auto & BE = *B.data();
  auto & CE = *C.data();
  for (std::size_t i=0; i<A.size(); ++i)
  {
    const auto & currMatrixEigen = *(A(i).data());
    assert( C.extent(0) == currMatrixEigen.rows() );
    assert( B.extent(0) == currMatrixEigen.cols() );
    CE.col(i) = beta * CE.col(i) + alpha * currMatrixEigen * BE.col(i);
  }
}

/***********************************
* special case A==B and op(A) = transpose
**********************************/
template <typename A_type, typename scalar_type, typename C_type>
::pressio::mpl::enable_if_t<
  ::pressio::ops::constraints::rank2_container_eigen_with_native_data_access<A_type>::value and
  ::pressio::ops::constraints::rank2_container_eigen_with_native_data_access<C_type>::value
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const A_type & A,
	const scalar_type beta,
	C_type & C)
{
  static_assert
    (containers::predicates::are_scalar_compatible<A_type, C_type>::value,
     "Types are not scalar compatible");

  const auto & AE = *A.data();
  auto & CE = *C.data();
  CE = beta * CE + alpha * AE.transpose() * AE;
}

template <typename C_type, typename A_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  ::pressio::ops::constraints::rank2_container_eigen_with_native_data_access<A_type>::value and
  ::pressio::ops::constraints::rank2_container_eigen_with_native_data_access<C_type>::value,
  C_type
  >
product(::pressio::transpose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const A_type & A)
{
  static_assert
    (containers::predicates::are_scalar_compatible<A_type, C_type>::value,
     "Types are not scalar compatible");

  constexpr auto zero = ::pressio::utils::constants<scalar_type>::zero();
  C_type C(A.extent(1), A.extent(1));
  product(modeA, modeB, alpha, A, A, zero, C);
  return C;
}

//-------------------------------------------
// C = beta * C + alpha*A*B
// specialize for when A = asDiagonalMatrix expression
//-------------------------------------------
template <typename T, typename B_type, typename scalar_type, typename C_type>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_vector_wrapper_eigen<T>::value and
  ::pressio::ops::constraints::rank2_container_eigen_with_native_data_access<B_type>::value and
  ::pressio::ops::constraints::rank2_container_eigen_with_native_data_access<C_type>::value
  >
product(::pressio::nontranspose modeA,
	::pressio::nontranspose modeB,
	const scalar_type alpha,
	const pressio::containers::expressions::AsDiagonalMatrixExpr<T> & A,
	const B_type & B,
	const scalar_type beta,
	C_type & C)
{
  static_assert
    (containers::predicates::are_scalar_compatible<T, B_type, C_type>::value,
     "Types are not scalar compatible");

  assert( C.extent(0) == A.extent(0) );
  assert( C.extent(1) == B.extent(1) );
  assert( A.extent(1) == B.extent(0) );
  const auto & AE = *(A.pressioObj()->data());
  const auto & BE = *B.data();
  auto & CE = *C.data();
  CE = beta*CE + alpha * (AE.asDiagonal() * BE);
}

}}//end namespace pressio::ops
#endif  // OPS_EIGEN_OPS_LEVEL3_HPP_
