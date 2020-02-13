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

#include "KokkosBlas2_gemv.hpp"

namespace pressio{ namespace containers{ namespace ops{

/*
 * multi_vector prod vector
 */

// begin namespace pressio::containers::ops::impl
namespace impl{

template <
  typename mvec_type, typename operand_t, typename res_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_tpetra<mvec_type>::value and
    containers::meta::is_vector_wrapper_tpetra<res_type>::value and
    !containers::meta::is_vector_wrapper_kokkos<operand_t>::value
    > * = nullptr
  >
void _product_tpetra_mv_sharedmem_vec(const mvec_type & mvA,
				      const operand_t & b,
				      res_type & C){

  //zero out result
  ::pressio::containers::ops::set_zero(C);

  // how many vectors are in mvA
  const auto numVecs = mvA.numVectors();
  // size of vecB
  assert(size_t(numVecs) == size_t(b.extent(0)));

  // my number of rows
  const auto myNrows = mvA.extentLocal(0);

  // get the wrapped trilinos tpetra multivector
  auto trilD = mvA.data();
  //  trilD->template sync<Kokkos::HostSpace>();
  auto mv2d = trilD->template getLocalView<Kokkos::HostSpace>();

  // get wrapped data for the result too
  auto C1 = C.data()->template getLocalView<Kokkos::HostSpace>();
  auto C2 = Kokkos::subview(C1, Kokkos::ALL(), 0);
  C.data()->template modify<Kokkos::HostSpace>();

  // loop
  for (size_t i=0; i<(size_t)myNrows; i++){
    for (size_t j=0; j<(size_t)numVecs; j++){
     C2[i] += mv2d(i,j) * b[j];
    }
  }
  using device_t = typename details::traits<res_type>::device_t;
  C.data()->template sync<device_t>();
}


// when the operand is a kokkos wrapper we use kokkos functionalities directly
template <
  typename mvec_type, typename operand_t, typename res_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_tpetra<mvec_type>::value and
    containers::meta::is_vector_wrapper_kokkos<operand_t>::value and
    containers::meta::is_vector_wrapper_tpetra<res_type>::value
    > * = nullptr
  >
void _product_tpetra_mv_sharedmem_vec(const mvec_type & mvA,
				      const operand_t & b,
				      res_type & C){

  // make sure the tpetra mv has same exe space of the kokkos vector wrapper
  using tpetra_mv_dev_t = typename ::pressio::containers::details::traits<mvec_type>::device_t;
  using kokkos_v_dev_t  = typename ::pressio::containers::details::traits<operand_t>::device_type;
  static_assert( std::is_same<tpetra_mv_dev_t, kokkos_v_dev_t>::value,
		 "product: tpetra MV and kokkos wrapper need to have same device type" );
  using dev_t  = tpetra_mv_dev_t;

  assert( mvA.numVectors() == b.data()->extent(0) );

  using sc_t = typename containers::details::traits<mvec_type>::scalar_t;
  constexpr auto zero = ::pressio::utils::constants::zero<sc_t>();
  constexpr auto one = ::pressio::utils::constants::one<sc_t>();
  const char ctA = 'N';

  const auto mvALocalView_d = mvA.data()->template getLocalView<dev_t>();
  // I need to do the following because Tpetra::Vector is implemented
  // as a special case of MultiVector so getLocalView returns a rank-2 view
  // so in order to get view with rank==1 I need to explicitly get the subview
  const auto mvCLocalView_drank2 = C.data()->template getLocalView<dev_t>();
  const auto mvCLocalView_drank1 = Kokkos::subview(mvCLocalView_drank2, Kokkos::ALL(), 0);
  KokkosBlas::gemv(&ctA, one, mvALocalView_d, *b.data(), zero, mvCLocalView_drank1);
}

}//end namespace pressio::containers::ops::impl




/* -------------------------------------------------------------------
 * specialize for tpetra mv operating on a sharedmem vector wrapper
 *-------------------------------------------------------------------*/
template <
  typename mvec_type,
  typename vec_type,
  typename res_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_tpetra<mvec_type>::value and
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value and
    containers::meta::is_vector_wrapper_tpetra<res_type>::value and
    (containers::meta::is_vector_wrapper_eigen<vec_type>::value or
     containers::meta::is_dense_vector_wrapper_teuchos<vec_type>::value or
     containers::meta::is_vector_wrapper_kokkos<vec_type>::value)
    > * = nullptr
  >
void product(const mvec_type & mvA, const vec_type & vecB, res_type & C)
{
  ::pressio::containers::ops::impl::_product_tpetra_mv_sharedmem_vec(mvA, vecB, C);
}

template <
  typename mvec_type,
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_tpetra<mvec_type>::value and
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value and
    (containers::meta::is_vector_wrapper_eigen<vec_type>::value or
     containers::meta::is_dense_vector_wrapper_teuchos<vec_type>::value or
     containers::meta::is_vector_wrapper_kokkos<vec_type>::value)
    > * = nullptr
 >
auto product(const mvec_type & mvA, const vec_type & vecB)
  -> containers::Vector<
  Tpetra::Vector<typename details::traits<mvec_type>::scalar_t,
                 typename details::traits<mvec_type>::local_ordinal_t,
                 typename details::traits<mvec_type>::global_ordinal_t,
                 typename details::traits<mvec_type>::node_t>
                 >
{

  // here, mvA is distrubted, but vecB is NOT.

  // the data map of the multivector
  auto rcpMap = mvA.data()->getMap();

  using mvec_traits = typename details::traits<mvec_type>;
  using sc_t = typename mvec_traits::scalar_t;
  using LO_t = typename mvec_traits::local_ordinal_t;
  using GO_t = typename mvec_traits::global_ordinal_t;
  using NO_t = typename mvec_traits::node_t;

  // result is an Tpetra Vector with same distribution of mvA
  using res_nat_t = Tpetra::Vector<sc_t, LO_t, GO_t, NO_t>;
  using res_t = containers::Vector<res_nat_t>;
  res_t c(rcpMap);
  product(mvA, vecB, c);
  return c;
}



/* -------------------------------------------------------------------
 * specialize for tpetra mv operating on an expression
 *-------------------------------------------------------------------*/
template <
  typename mvec_type,
  typename expr_type,
  typename res_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_tpetra<mvec_type>::value and
    ::pressio::containers::meta::is_expression<expr_type>::value
    > * = nullptr
  >
void product(const mvec_type & mvA, const expr_type & b, res_type & C)
{
  ::pressio::containers::ops::impl::_product_tpetra_mv_sharedmem_vec(mvA, b, C);
}

template <
  typename mvec_type,
  typename expr_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_tpetra<mvec_type>::value and
    ::pressio::containers::meta::is_expression<expr_type>::value
    > * = nullptr
 >
auto product(const mvec_type & mvA, const expr_type & b)
  -> containers::Vector<
  Tpetra::Vector<typename details::traits<mvec_type>::scalar_t,
                 typename details::traits<mvec_type>::local_ordinal_t,
                 typename details::traits<mvec_type>::global_ordinal_t,
                 typename details::traits<mvec_type>::node_t>
                 >
{
  // the data map of the multivector
  auto rcpMap = mvA.data()->getMap();

  using mvec_traits = typename details::traits<mvec_type>;
  using sc_t = typename mvec_traits::scalar_t;
  using LO_t = typename mvec_traits::local_ordinal_t;
  using GO_t = typename mvec_traits::global_ordinal_t;
  using NO_t = typename mvec_traits::node_t;

  // result is an Tpetra Vector with same distribution of mvA
  using res_nat_t = Tpetra::Vector<sc_t, LO_t, GO_t, NO_t>;
  using res_t = containers::Vector<res_nat_t>;
  res_t c(rcpMap);
  product(mvA, b, c);
  return c;
}


}}}//end namespace pressio::containers::ops
#endif
#endif
