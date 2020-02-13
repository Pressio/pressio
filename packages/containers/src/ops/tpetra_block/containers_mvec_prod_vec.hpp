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
#ifndef CONTAINERS_SRC_OPS_TPETRA_BLOCK_MULTI_VECTOR_PROD_VECTOR_HPP_
#define CONTAINERS_SRC_OPS_TPETRA_BLOCK_MULTI_VECTOR_PROD_VECTOR_HPP_

namespace pressio{ namespace containers{ namespace ops{

/*
 * tpetra block multi_vector prod shared mem vector
 */

/* -------------------------------------------------------------------
 * specialize for tpetra block mv product Kokkos wrapper
 *-------------------------------------------------------------------*/
template <
  typename mvec_type,
  typename vec_type,
  typename res_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_tpetra_block<mvec_type>::value and
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value and
    containers::meta::is_vector_wrapper_kokkos<vec_type>::value and
    containers::meta::is_vector_wrapper_tpetra_block<res_type>::value
    > * = nullptr
  >
void product(const mvec_type & mvA, const vec_type & vecB, res_type & C)
{
  using sc_t = typename containers::details::traits<mvec_type>::scalar_t;

  // make sure the tpetra mv has same exe space of the kokkos vector wrapper
  using tpetra_mv_dev_t = typename ::pressio::containers::details::traits<mvec_type>::device_t;
  using kokkos_v_dev_t  = typename ::pressio::containers::details::traits<vec_type>::device_type;
  static_assert( std::is_same<tpetra_mv_dev_t, kokkos_v_dev_t>::value,
		 "product: tpetra MV and kokkos wrapper need to have same device type" );

  assert( mvA.globalNumVectors() == vecB.size() );

  constexpr auto zero = ::pressio::utils::constants::zero<sc_t>();
  constexpr auto one = ::pressio::utils::constants::one<sc_t>();
  const char ctA = 'N';

  // the the underlying tpetra multivector
  const auto mvView = mvA.data()->getMultiVectorView();
  // get a local view
  const auto mvALocalView_d = mvView.getLocalViewDevice();
  // I need to do the following because Tpetra::Vector is implemented
  // as a special case of MultiVector so getLocalView returns a rank-2 view
  // so in order to get view with rank==1 I need to explicitly get the subview
  const auto CView = C.data()->getVectorView();
  const auto CLocalView_drank2 = CView.getLocalViewDevice();
  const auto CLocalView_drank1 = Kokkos::subview(CLocalView_drank2, Kokkos::ALL(), 0);
  KokkosBlas::gemv(&ctA, one, mvALocalView_d, *vecB.data(), zero, CLocalView_drank1);
}

template<
  typename mvec_type,
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_tpetra_block<mvec_type>::value and
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value and
    containers::meta::is_vector_wrapper_kokkos<vec_type>::value
    > * = nullptr
  >
auto product(const mvec_type & mvA, const vec_type & vecB)
  -> containers::Vector<
    Tpetra::Experimental::BlockVector<
      typename details::traits<mvec_type>::scalar_t,
      typename details::traits<mvec_type>::local_ordinal_t,
      typename details::traits<mvec_type>::global_ordinal_t,
      typename details::traits<mvec_type>::node_t
      >
    >
{
  // the data map of the multivector
  const auto rcpMap = mvA.getRCPDataMap();
  // the block size
  const auto mvABlockSize = mvA.getBlockSize();

  using mvec_traits = typename details::traits<mvec_type>;
  using sc_t = typename mvec_traits::scalar_t;
  using LO_t = typename mvec_traits::local_ordinal_t;
  using GO_t = typename mvec_traits::global_ordinal_t;
  using NO_t = typename mvec_traits::node_t;

  // result is an Tpetra Vector with same distribution of mvA
  using res_nat_t = Tpetra::Experimental::BlockVector<sc_t, LO_t, GO_t, NO_t>;
  using res_t = containers::Vector<res_nat_t>;
  res_t c( res_nat_t(*rcpMap, mvABlockSize) );
  product(mvA, vecB, c);
  return c;
}



/* -------------------------------------------------------------------
 * specialize for tpetra block mv product eigen wrapper
 *-------------------------------------------------------------------*/
template <
  typename mvec_type,
  typename vec_type,
  typename res_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_tpetra_block<mvec_type>::value and
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value and
    containers::meta::is_vector_wrapper_tpetra_block<res_type>::value and
    containers::meta::is_vector_wrapper_eigen<vec_type>::value
    > * = nullptr
  >
void product(const mvec_type & mvA, const vec_type & vecB, res_type & C)
{
  /* computes: C = mvA*vecB */

  //zero out result
  ::pressio::containers::ops::set_zero(C);

  // how many vectors are in mvA
  const auto numVecs = mvA.globalNumVectors();

  // size of vecB
  const size_t vecBLen = vecB.size();
  if (vecBLen != size_t(numVecs))
    assert(size_t(numVecs) == vecBLen);

  // get the wrapped trilinos tpetra multivector
  const auto mvA_mvv = mvA.data()->getMultiVectorView();
  const auto mvA_hv = mvA_mvv.template getLocalView<Kokkos::HostSpace>();
  // my number of rows
  const auto myNrows = mvA_mvv.getLocalLength();

  // the result is a block tpetra vector, get the regular tpetra vector
  auto C_vv = C.data()->getVectorView();
  auto C_hv = C_vv.template getLocalView<Kokkos::HostSpace>();
  C_vv.template modify<Kokkos::HostSpace>();

  for (std::size_t i=0; i<(std::size_t)myNrows; i++){
    for (std::size_t j=0; j<(std::size_t)numVecs; j++){
      // we use C_hv(i,0) because C is Tpetra::Vector, which is
      // actually a Tpetra::MultiVector with one column, so C_hv
      // is a kokkos::View<scalar**,...> so we need to index the zero column
      C_hv(i,0) += mvA_hv(i,j) * vecB[j];
    }
  }
  using device_t = typename details::traits<res_type>::device_t;
  C.data()->template sync<device_t>();
}

template<
  typename mvec_type,
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_tpetra_block<mvec_type>::value and
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value and
    containers::meta::is_vector_wrapper_eigen<vec_type>::value
    > * = nullptr
  >
auto product(const mvec_type & mvA, const vec_type & vecB)
  -> containers::Vector<
    Tpetra::Experimental::BlockVector<
      typename details::traits<mvec_type>::scalar_t,
      typename details::traits<mvec_type>::local_ordinal_t,
      typename details::traits<mvec_type>::global_ordinal_t,
      typename details::traits<mvec_type>::node_t
      >
    >
{
  // the data map of the multivector
  const auto rcpMap = mvA.getRCPDataMap();
  // the block size
  const auto mvABlockSize = mvA.getBlockSize();

  using mvec_traits = typename details::traits<mvec_type>;
  using sc_t = typename mvec_traits::scalar_t;
  using LO_t = typename mvec_traits::local_ordinal_t;
  using GO_t = typename mvec_traits::global_ordinal_t;
  using NO_t = typename mvec_traits::node_t;

  // result is an Tpetra Vector with same distribution of mvA
  using res_nat_t = Tpetra::Vector<sc_t, LO_t, GO_t, NO_t>;
  using res_t = containers::Vector<res_nat_t>;
  res_t c( res_nat_t(*rcpMap, mvABlockSize) );
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
    containers::meta::is_multi_vector_wrapper_tpetra_block<mvec_type>::value and
    ::pressio::containers::meta::is_expression<expr_type>::value
    > * = nullptr
  >
void product(const mvec_type & mvA, const expr_type & b, res_type & C)
{
  throw std::runtime_error("Warning, container::ops:product operation between tpetra block and expression not yet supported"); 
  //::pressio::containers::ops::impl::_product_tpetra_mv_sharedmem_vec(mvA, b, C);
}

template <
  typename mvec_type,
  typename expr_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_tpetra_block<mvec_type>::value and
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
  throw std::runtime_error("Warning, container::ops::product operation between tpetra block and expression not yet supported");
  // the data map of the multivector
  /* 
  auto rcpMap = mvA.getRCPDataMap();

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
  */
}


}}}//end namespace pressio::containers::ops
#endif
#endif
