/*
//@HEADER
// ************************************************************************
//
// containers_mvec_prod_vec.hpp
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

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#ifndef CONTAINERS_SRC_OPS_TPETRA_MULTI_VECTOR_PROD_VECTOR_HPP_
#define CONTAINERS_SRC_OPS_TPETRA_MULTI_VECTOR_PROD_VECTOR_HPP_

#include "Tpetra_idot.hpp"
#include "KokkosBlas2_gemv.hpp"

namespace pressio{ namespace containers{ namespace ops{


/*
 * multi_vector prod vector
 *
 * y = beta * y + alpha*op(A)*x
 *
*/

// begin namespace pressio::containers::ops::impl
namespace impl{

template <
  typename A_type, typename x_type, typename y_type, typename scalar_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_tpetra<A_type>::value and
    containers::meta::is_vector_wrapper_tpetra<y_type>::value
    > * = nullptr
  >
void _product_tpetra_mv_sharedmem_vec(const scalar_type alpha,
				      const A_type & A,
				      const x_type & x,
				      const scalar_type beta,
				      y_type & y)
{
  static_assert(containers::meta::are_scalar_compatible<A_type, x_type, y_type>::value,
		"Types are not scalar compatible");

  // how many vectors are in A
  const auto numVecs = A.numVectors();
  // size of vecB
  assert(size_t(numVecs) == size_t(x.extent(0)));

  // my number of rows
  const auto myNrows = A.extentLocal(0);

  // get the wrapped trilinos tpetra multivector
  auto trilD = A.data();
  //  trilD->template sync<Kokkos::HostSpace>();
  auto mv2d = trilD->template getLocalView<Kokkos::HostSpace>();

  // get wrapped data for the result too
  auto y1 = y.data()->template getLocalView<Kokkos::HostSpace>();
  auto y2 = Kokkos::subview(y1, Kokkos::ALL(), 0);
  y.data()->template modify<Kokkos::HostSpace>();

  // loop
  for (size_t i=0; i<(size_t)myNrows; i++){
    y2[i] = beta*y2[i];
    for (size_t j=0; j<(size_t)numVecs; j++){
      y2[i] += alpha * mv2d(i,j) * x[j];
    }
  }
  using device_t = typename details::traits<y_type>::device_t;
  y.data()->template sync<device_t>();
}


// when the operand is a kokkos wrapper we use kokkos functionalities directly
template <
  typename A_type, typename x_type, typename y_type, typename scalar_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_tpetra<A_type>::value and
    containers::meta::is_vector_wrapper_tpetra<y_type>::value and
    containers::meta::is_vector_wrapper_kokkos<x_type>::value
    > * = nullptr
  >
void _product_tpetra_mv_sharedmem_vec_kokkos(const scalar_type alpha,
					     const A_type & A,
					     const x_type & x,
					     const scalar_type beta,
					     y_type & y)
{
  static_assert(containers::meta::are_scalar_compatible<A_type, x_type, y_type>::value,
    "Types are not scalar compatible");

  // make sure the tpetra mv has same exe space of the kokkos vector wrapper
  using tpetra_mv_dev_t = typename ::pressio::containers::details::traits<A_type>::device_t;
  using kokkos_v_dev_t  = typename ::pressio::containers::details::traits<x_type>::device_type;
  static_assert( std::is_same<tpetra_mv_dev_t, kokkos_v_dev_t>::value,
		 "product: tpetra MV and kokkos wrapper need to have same device type" );
  using dev_t  = tpetra_mv_dev_t;

  assert( A.numVectors() == x.data()->extent(0) );

  using sc_t = typename containers::details::traits<A_type>::scalar_t;
  // constexpr auto zero = ::pressio::utils::constants::zero<sc_t>();
  // constexpr auto one = ::pressio::utils::constants::one<sc_t>();
  const char ctA = 'N';

  const auto ALocalView_d = A.data()->template getLocalView<dev_t>();
  // I need to do the following because Tpetra::Vector is implemented
  // as a special case of MultiVector so getLocalView returns a rank-2 view
  // so in order to get view with rank==1 I need to explicitly get the subview
  const auto mvCLocalView_drank2 = y.data()->template getLocalView<dev_t>();
  const auto mvCLocalView_drank1 = Kokkos::subview(mvCLocalView_drank2, Kokkos::ALL(), 0);
  KokkosBlas::gemv(&ctA, alpha, ALocalView_d, *x.data(), beta, mvCLocalView_drank1);
}

}//end namespace pressio::containers::ops::impl



/* -------------------------------------------------------------------
 * x is a sharedmem vector wrapper
 *-------------------------------------------------------------------*/
template < typename A_type, typename x_type, typename scalar_type, typename y_type>
::pressio::mpl::enable_if_t<
  containers::meta::is_multi_vector_wrapper_tpetra<A_type>::value and
  containers::meta::is_vector_wrapper_tpetra<y_type>::value and
  (containers::meta::is_vector_wrapper_eigen<x_type>::value or
   containers::meta::is_dense_vector_wrapper_teuchos<x_type>::value)
  >
product(::pressio::nontranspose mode,
	const scalar_type alpha,
	const A_type & A,
	const x_type & x,
	const scalar_type beta,
	y_type & y)
{
  ::pressio::containers::ops::impl::_product_tpetra_mv_sharedmem_vec(alpha, A, x, beta, y);
}

template < typename A_type, typename x_type, typename scalar_type, typename y_type>
::pressio::mpl::enable_if_t<
  containers::meta::is_multi_vector_wrapper_tpetra<A_type>::value and
  containers::meta::is_vector_wrapper_tpetra<y_type>::value and
  containers::meta::is_vector_wrapper_kokkos<x_type>::value
  >
product(::pressio::nontranspose mode,
	const scalar_type alpha,
	const A_type & A,
	const x_type & x,
	const scalar_type beta,
	y_type & y)
{
  ::pressio::containers::ops::impl::_product_tpetra_mv_sharedmem_vec_kokkos(alpha, A, x, beta, y);
}



/* -------------------------------------------------------------------
 * x is a distributed Tpetra vector wrapper
 *-------------------------------------------------------------------*/

// y = wrapper of Kokkos vector
template <typename A_type, typename x_type, typename y_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  containers::meta::is_multi_vector_wrapper_tpetra<A_type>::value and
  containers::meta::is_vector_wrapper_tpetra<x_type>::value and
  containers::meta::is_vector_wrapper_kokkos<y_type>::value
  >
product(::pressio::transpose mode,
	const scalar_type alpha,
	const A_type & A,
	const x_type & x,
	const scalar_type beta,
	y_type & y)
{
  static_assert(containers::meta::are_scalar_compatible<A_type, x_type, y_type>::value,
		"Tpetra MV dot V: operands do not have matching scalar type");

  static_assert(std::is_same<
		typename containers::details::traits<A_type>::device_t,
		typename containers::details::traits<x_type>::device_t>::value,
		"Tpetra MV dot V: operands do not have the same device type");

  static_assert(std::is_same<
		typename containers::details::traits<x_type>::device_t,
		typename containers::details::traits<y_type>::device_t>::value,
		"Tpetra MV dot V: V and result do not have the same device type");

  auto request = Tpetra::idot( *y.data(), *A.data(), *x.data());
  request->wait();
}


// y = scalar *, passed in
template <typename A_type, typename x_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  containers::meta::is_multi_vector_wrapper_tpetra<A_type>::value and
  containers::meta::is_vector_wrapper_tpetra<x_type>::value
  >
product(::pressio::transpose mode,
	const scalar_type alpha,
	const A_type & A,
	const x_type & x,
	const scalar_type beta,
	scalar_type * y)
{
  static_assert(containers::meta::are_scalar_compatible<A_type, x_type>::value,
		"Types are not scalar compatible");

  const auto numVecs = A.numVectors();
  for (auto i=0; i<numVecs; i++){
    // colI is a Teuchos::RCP<Vector<...>>
    const auto colI = A.data()->getVector(i);
    y[i] = beta * y[i] + alpha * colI->dot(*x.data());
  }
}

// y = Eigen vector
template <typename A_type, typename x_type, typename y_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  containers::meta::is_multi_vector_wrapper_tpetra<A_type>::value and
  containers::meta::is_vector_wrapper_tpetra<x_type>::value and
  containers::meta::is_vector_wrapper_eigen<y_type>::value
  >
product(::pressio::transpose mode,
	const scalar_type alpha,
	const A_type & A,
	const x_type & x,
	const scalar_type beta,
	y_type & y)
{
  ///computes dot product of each vector in A
  ///with vecB storing each value in result

  /* Apparently, trilinos does not support this...
     the native dot() method of multivectors is only for
     dot product of two multivectors with same # of columns.
     So we have to extract each column vector
     from A and do dot product one a time*/

  static_assert(containers::meta::are_scalar_compatible<A_type, x_type, y_type>::value,
		"Types are not scalar compatible");

  const auto numVecs = A.numVectors();
  // check the result has right size
  assert( y.extent(0) == numVecs );
  product(mode, alpha, A, x, beta, y.data()->data());
}

}}}//end namespace pressio::containers::ops
#endif
#endif





// /* -------------------------------------------------------------------
//  * specialize for tpetra mv operating on an expression
//  *-------------------------------------------------------------------*/
// template <
//   typename mvec_type,
//   typename expr_type,
//   typename res_type,
//   ::pressio::mpl::enable_if_t<
//     containers::meta::is_multi_vector_wrapper_tpetra<mvec_type>::value and
//     ::pressio::containers::meta::is_expression<expr_type>::value
//     > * = nullptr
//   >
// void product(const mvec_type & mvA, const expr_type & b, res_type & C)
// {
//   ::pressio::containers::ops::impl::_product_tpetra_mv_sharedmem_vec(mvA, b, C);
// }

// template <
//   typename mvec_type,
//   typename expr_type,
//   ::pressio::mpl::enable_if_t<
//     containers::meta::is_multi_vector_wrapper_tpetra<mvec_type>::value and
//     ::pressio::containers::meta::is_expression<expr_type>::value
//     > * = nullptr
//  >
// auto product(const mvec_type & mvA, const expr_type & b)
//   -> containers::Vector<
//   Tpetra::Vector<typename details::traits<mvec_type>::scalar_t,
//                  typename details::traits<mvec_type>::local_ordinal_t,
//                  typename details::traits<mvec_type>::global_ordinal_t,
//                  typename details::traits<mvec_type>::node_t>
//                  >
// {
//   // the data map of the multivector
//   auto rcpMap = mvA.data()->getMap();

//   using mvec_traits = typename details::traits<mvec_type>;
//   using sc_t = typename mvec_traits::scalar_t;
//   using LO_t = typename mvec_traits::local_ordinal_t;
//   using GO_t = typename mvec_traits::global_ordinal_t;
//   using NO_t = typename mvec_traits::node_t;

//   // result is an Tpetra Vector with same distribution of mvA
//   using res_nat_t = Tpetra::Vector<sc_t, LO_t, GO_t, NO_t>;
//   using res_t = containers::Vector<res_nat_t>;
//   res_t c(rcpMap);
//   product(mvA, b, c);
//   return c;
// }




// //--------------------------------------
// // compute y = y + mvA^T vecB
// // y = eigen expression
// //--------------------------------------
// template <
//   typename mvec_type,
//   typename vec_type,
//   typename expr_type,
//   ::pressio::mpl::enable_if_t<
//     containers::meta::is_multi_vector_wrapper_tpetra<mvec_type>::value and
//     containers::meta::is_vector_wrapper_tpetra<vec_type>::value and
//     containers::meta::is_expression<expr_type>::value and
//     ::pressio::containers::meta::is_vector_wrapper_eigen<
//       typename ::pressio::containers::details::traits<expr_type>::data_t
//       >::value
//     > * = nullptr
//   >
// void updateWithDot(const mvec_type & mvA, const vec_type & vecB, expr_type & result)
// {
//   static_assert(containers::meta::are_scalar_compatible<mvec_type, vec_type, expr_type>::value,
//     "Types are not scalar compatible");

//   // check the result has right size
//   const auto numVecs = mvA.numVectors();
//   assert( result.extent(0) == numVecs );

//   for (auto i=0; i<numVecs; i++){
//     // colI is a Teuchos::RCP<Vector<...>>
//     auto colI = mvA.data()->getVector(i);
//     result[i] += colI->dot(*vecB.data());
//   }
// }



// //--------------------------------------------
// // c = teuchos serial dense vector, passed in
// //--------------------------------------------
// template <
//   typename mvec_type,
//   typename vec_type,
//   typename result_vec_type,
//   ::pressio::mpl::enable_if_t<
//     containers::meta::is_multi_vector_wrapper_tpetra<mvec_type>::value and
//     containers::meta::is_vector_wrapper_tpetra<vec_type>::value and
//     containers::meta::is_dense_vector_wrapper_teuchos<result_vec_type>::value and
//     containers::details::traits<result_vec_type>::is_dynamic
//     > * = nullptr
//   >
// void dot(const mvec_type & mvA,
// 	 const vec_type & vecB,
// 	 result_vec_type & result)
// {
//   static_assert(containers::meta::are_scalar_compatible<mvec_type, vec_type, result_vec_type>::value,
//     "Types are not scalar compatible");

//   const auto numVecs = mvA.numVectors();
//   if ( result.extent(0) != numVecs )
//     result.data()->resize(numVecs);
//   dot(mvA, vecB, result.data()->values());
// }
