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

#ifndef OPS_EPETRA_OPS_LEVEL2_HPP_
#define OPS_EPETRA_OPS_LEVEL2_HPP_

namespace pressio{ namespace ops{

/*
 * y = beta * y + alpha*op(A)*x
*/

// begin namespace pressio::ops::impl
namespace impl{

template <typename A_type, typename x_type, typename scalar_type>
void _product_epetra_mv_sharedmem_vec(const scalar_type alpha,
				      const A_type & A,
				      const x_type & x,
				      const scalar_type beta,
				      Epetra_Vector & y)
{
  constexpr auto zero = pressio::utils::Constants<scalar_type>::zero();
  ::pressio::ops::scale(y, beta);
  if (alpha == zero) {
    return;
  }

  const int numVecs = A.NumVectors();
  for (int i=0; i< A.MyLength(); i++)
  {
    for (int j=0; j< (int)numVecs; j++){
      y[i] += alpha * A[j][i] * x(j);
    }
  }
}
}//end namespace pressio::ops::impl

// -------------------------------
// y = beta * y + alpha*A*x
//
// y : epetra vector
// A : epetra MultiVector
// x : Teuchos Vector
// -------------------------------
template <
  typename alpha_t,
  typename A_type,
  typename x_type,
  typename y_type,
  typename beta_t
  >
std::enable_if_t<
  // level2 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<x_type>::rank == 1
  && ::pressio::Traits<y_type>::rank == 1
  // TPL/container specific
  && ::pressio::is_multi_vector_epetra<A_type>::value
  && ::pressio::is_dense_vector_teuchos<x_type>::value
  && ::pressio::is_vector_epetra<y_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, x_type, y_type>::value
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  >
product(::pressio::nontranspose /*unused*/,
	const alpha_t alpha,
	const A_type & A,
	const x_type & x,
	const beta_t beta,
	y_type & y)
{
  ::pressio::ops::impl::_product_epetra_mv_sharedmem_vec(alpha, A, x, beta, y);
}

// -------------------------------
// y = beta * y + alpha*A^T*x
//
// y : Teuchos vector
// A : epetra MultiVector
// x : epetra Vector
// -------------------------------
template <
  typename alpha_t,
  typename A_type,
  typename x_type,
  typename y_type,
  typename beta_t
  >
std::enable_if_t<
  // level2 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<x_type>::rank == 1
  && ::pressio::Traits<y_type>::rank == 1
  // TPL/container specific
  && ::pressio::is_multi_vector_epetra<A_type>::value
  && ::pressio::is_vector_epetra<x_type>::value
  && ::pressio::is_dense_vector_teuchos<y_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, x_type, y_type>::value
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  >
product(::pressio::transpose /*unused*/,
  const alpha_t alpha,
  const A_type & A,
  const x_type & x,
  const beta_t beta,
  y_type & y)
{

  const int numVecs = A.NumVectors();
  assert( (std::size_t)y.length() == (std::size_t)numVecs );

  using sc_t = typename ::pressio::Traits<A_type>::scalar_type;
  const auto zero = ::pressio::utils::Constants<sc_t>::zero();
  auto tmp = zero;
  for (int i=0; i<numVecs; i++)
  {
    y(i) = (beta == zero) ? zero : beta * y(i);
    if (!(alpha == zero)) {
      A(i)->Dot(x, &tmp);
      y(i) += alpha * tmp;
    }
  }
}

#ifdef PRESSIO_ENABLE_TPL_EIGEN
// -------------------------------
// y = beta * y + alpha*A*x
//
// y : epetra vector
// A : epetra MultiVector
// x : eigen Vector
// -------------------------------
template <
  typename alpha_t,
  typename A_type,
  typename x_type,
  typename y_type,
  typename beta_t
  >
std::enable_if_t<
  // level2 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<x_type>::rank == 1
  && ::pressio::Traits<y_type>::rank == 1
  // TPL/container specific
  && ::pressio::is_multi_vector_epetra<A_type>::value
  && (::pressio::is_vector_eigen<x_type>::value
   || ::pressio::is_expression_acting_on_eigen<x_type>::value)
  && ::pressio::is_vector_epetra<y_type>::value
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<A_type, x_type, y_type>::value
  && std::is_convertible<alpha_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && std::is_convertible<beta_t, typename ::pressio::Traits<A_type>::scalar_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  >
product(::pressio::nontranspose /*unused*/,
  const alpha_t alpha,
  const A_type & A,
  const x_type & x,
  const beta_t beta,
  y_type & y)
{
  ::pressio::ops::impl::_product_epetra_mv_sharedmem_vec(alpha, A, x, beta, y);
}

// -------------------------------
// y = beta * y + alpha*A^T*x
//
// y : eigen vector
// A : epetra MultiVector
// x : epetra Vector
// -------------------------------
template <
  typename alpha_t,
  typename A_type,
  typename x_type,
  typename y_type,
  typename beta_t
  >
std::enable_if_t<
  // level2 common constraints
     ::pressio::Traits<A_type>::rank == 2
  && ::pressio::Traits<x_type>::rank == 1
  && ::pressio::Traits<y_type>::rank == 1
  // TPL/container specific
  && ::pressio::is_multi_vector_epetra<A_type>::value
  && ::pressio::is_vector_epetra<x_type>::value
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
  const alpha_t alpha,
  const A_type & A,
  const x_type & x,
  const beta_t beta,
  y_type & y)
{

  const int numVecs = A.NumVectors();
  assert( (std::size_t)::pressio::ops::extent(y, 0) == (std::size_t)numVecs );

  using sc_t = typename ::pressio::Traits<A_type>::scalar_type;
  const auto zero = ::pressio::utils::Constants<sc_t>::zero();
  auto tmp = zero;
  for (int i=0; i<numVecs; i++)
  {
    y(i) = (beta == zero) ? zero : beta * y(i);
    if (!(alpha == zero)) {
      A(i)->Dot(x, &tmp);
      y(i) += alpha * tmp;
    }
  }
}
#endif


// #ifdef PRESSIO_ENABLE_TPL_EIGEN
// /* -------------------------------------------------------------------
//  * y = beta * y + alpha*A*x
//  * x is a eigen
//  * A = Eigen::MultiVector
//  * y = Epetra Vector
//  * this covers the case where the matrix A acts locally to x
//    while y is distributed, so A*x only fills the corresponding part of y
//  *-------------------------------------------------------------------*/
// template < typename A_type, typename x_type, typename scalar_type>
// std::enable_if_t<
//   ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<A_type>::value and
//   ::pressio::containers::predicates::is_vector_wrapper_eigen<x_type>::value
//   >
// product(::pressio::nontranspose mode,
// 	const scalar_type alpha,
// 	const A_type & A,
// 	const x_type & x,
// 	const scalar_type beta,
// 	::pressio::containers::Vector<Epetra_Vector> & y)
// {
//   static_assert
//     (std::is_same<
//      scalar_type, typename ::pressio::containers::details::traits<x_type>::scalar_t>::value,
//      "Scalar compatibility broken");

//   using eig_mapping_t = Eigen::Map< Eigen::Matrix<scalar_type, -1, -1> >;
//   eig_mapping_t yMapped(y.data()->Values(), y.extentLocal(0), 1);
//   const auto & AE = *A.data();
//   const auto & xE = *x.data();
//   yMapped = beta * yMapped + alpha * AE * xE;
// }
// #endif

// #ifdef PRESSIO_ENABLE_TPL_EIGEN
// /* -------------------------------------------------------------------
//  * y = beta * y + alpha*A^T*x
//  * x is a distributed Epetra vector wrapper
//  * A is an Eigen MV
//  * y is vector wrapper Eigen
//  *-------------------------------------------------------------------*/
// template <typename A_type, typename y_type, typename scalar_type>
// std::enable_if_t<
//   ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<A_type>::value and
//   ::pressio::containers::predicates::is_vector_wrapper_eigen<y_type>::value
//   >
// product(::pressio::transpose mode,
// 	const scalar_type alpha,
// 	const A_type & A,
// 	const ::pressio::containers::Vector<Epetra_Vector> & x,
// 	const scalar_type beta,
// 	y_type & y)
// {
//   using eig_mapping_t = const Eigen::Map< const Eigen::Matrix<scalar_type, -1, -1> >;
//   eig_mapping_t xMapped(x.data()->Values(), x.extentLocal(0), 1);
//   const auto & AE = *A.data();
//   auto & yE = *y.data();
//   yE = beta * yE + alpha * AE.transpose() * xMapped;
// }
// #endif

}}//end namespace pressio::ops
#endif  // OPS_EPETRA_OPS_LEVEL2_HPP_
