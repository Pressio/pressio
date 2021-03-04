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
 * multi_vector prod vector
 *
 * y = beta * y + alpha*op(A)*x
 *
*/

// begin namespace pressio::ops::impl
namespace impl{

template <typename A_type, typename x_type, typename scalar_type>
void _product_epetra_mv_sharedmem_vec(const scalar_type alpha,
				      const A_type & A,
				      const x_type & x,
				      const scalar_type beta,
				      containers::Vector<Epetra_Vector> & y)
{
  // how many vectors are in A
  const auto numVecs = A.numVectorsGlobal();

  // using go_t = typename ::pressio::containers::details::traits<A_type>::global_ordinal_t;
  // using lo_t = typename ::pressio::containers::details::traits<A_type>::local_ordinal_t;
  for (std::size_t i=0; i<(std::size_t)A.extentLocal(0); i++){
    //rowSum = zero;
    y(i) = beta*y(i);
    for (std::size_t j=0; j<(std::size_t)numVecs; j++){
      y(i) += alpha * A(i,j) * x(j);
    }
  }
}
}//end namespace pressio::ops::impl

/* -------------------------------------------------------------------
 * x is a sharedmem vector wrapper
 *-------------------------------------------------------------------*/
template < typename A_type, typename x_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_multi_vector_wrapper_epetra<A_type>::value and
  ::pressio::ops::constraints::sharedmem_host_subscriptable_rank1_container<x_type>::value
  >
product(::pressio::nontranspose mode,
	const scalar_type alpha,
	const A_type & A,
	const x_type & x,
	const scalar_type beta,
	::pressio::containers::Vector<Epetra_Vector> & y)
{
  static_assert
    (containers::predicates::are_scalar_compatible<A_type, x_type>::value,
     "Types are not scalar compatible");
  static_assert
    (mpl::is_same<
     scalar_type, typename ::pressio::containers::details::traits<x_type>::scalar_t>::value,
     "Scalar compatibility broken");

  ::pressio::ops::impl::_product_epetra_mv_sharedmem_vec(alpha, A, x, beta, y);
}


#ifdef PRESSIO_ENABLE_TPL_EIGEN
/* -------------------------------------------------------------------
 * y = beta * y + alpha*A*x
 * x is a sharedmem vector wrapper eigen
 * A = Eigen::MultiVector
 * y = Epetra Vector
 * this covers the case where the matrix A acts locally to x
   while y is distributed, so A*x only fills the corresponding part of y
 *-------------------------------------------------------------------*/
template < typename A_type, typename x_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<A_type>::value and
  ::pressio::containers::predicates::is_vector_wrapper_eigen<x_type>::value
  >
product(::pressio::nontranspose mode,
	const scalar_type alpha,
	const A_type & A,
	const x_type & x,
	const scalar_type beta,
	::pressio::containers::Vector<Epetra_Vector> & y)
{
  static_assert
    (containers::predicates::are_scalar_compatible<A_type, x_type>::value,
     "Types are not scalar compatible");
  static_assert
    (mpl::is_same<
     scalar_type, typename ::pressio::containers::details::traits<x_type>::scalar_t>::value,
     "Scalar compatibility broken");

  using eig_mapping_t = Eigen::Map< Eigen::Matrix<scalar_type, -1, -1> >;
  eig_mapping_t yMapped(y.data()->Values(), y.extentLocal(0), 1);
  const auto & AE = *A.data();
  const auto & xE = *x.data();
  yMapped = beta * yMapped + alpha * AE * xE;
}
#endif

/* -------------------------------------------------------------------
 * x is a distributed Epetra vector wrapper
 *-------------------------------------------------------------------*/
template <typename A_type, typename y_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_multi_vector_wrapper_epetra<A_type>::value and
  ::pressio::ops::constraints::sharedmem_host_subscriptable_rank1_container<y_type>::value
  >
product(::pressio::transpose mode,
	const scalar_type alpha,
	const A_type & A,
	const ::pressio::containers::Vector<Epetra_Vector> & x,
	const scalar_type beta,
	y_type & y)
{
  static_assert
    (::pressio::containers::predicates::are_scalar_compatible<A_type, y_type>::value,
     "Types are not scalar compatible");

  // using ord_t = typename ::pressio::containers::details::traits<A_type>::global_ordinal_t;
  const auto numVecs = A.numVectorsGlobal();
  assert( (std::size_t)y.extent(0) == (std::size_t)numVecs );

  auto * mvNatData = A.data();
  const auto * vecNatData = x.data();
  auto tmp = ::pressio::utils::constants<scalar_type>::zero();
  for (std::size_t i=0; i<(std::size_t)numVecs; i++)
  {
    (*mvNatData)(i)->Dot(*vecNatData, &tmp);
    y(i) = beta * y(i) + alpha * tmp;
  }
}

#ifdef PRESSIO_ENABLE_TPL_EIGEN
/* -------------------------------------------------------------------
 * y = beta * y + alpha*A^T*x
 * x is a distributed Epetra vector wrapper
 * A is an Eigen MV
 * y is vector wrapper Eigen
 *-------------------------------------------------------------------*/
template <typename A_type, typename y_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<A_type>::value and
  ::pressio::containers::predicates::is_vector_wrapper_eigen<y_type>::value
  >
product(::pressio::transpose mode,
	const scalar_type alpha,
	const A_type & A,
	const ::pressio::containers::Vector<Epetra_Vector> & x,
	const scalar_type beta,
	y_type & y)
{
  static_assert
    (::pressio::containers::predicates::are_scalar_compatible<A_type, y_type>::value,
     "Types are not scalar compatible");

  using eig_mapping_t = const Eigen::Map< const Eigen::Matrix<scalar_type, -1, -1> >;
  eig_mapping_t xMapped(x.data()->Values(), x.extentLocal(0), 1);
  const auto & AE = *A.data();
  auto & yE = *y.data();
  yE = beta * yE + alpha * AE.transpose() * xMapped;
}
#endif

}}//end namespace pressio::ops
#endif  // OPS_EPETRA_OPS_LEVEL2_HPP_
