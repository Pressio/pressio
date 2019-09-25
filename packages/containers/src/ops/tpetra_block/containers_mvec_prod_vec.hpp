
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_SRC_OPS_TPETRA_BLOCK_MULTI_VECTOR_PROD_VECTOR_HPP_
#define CONTAINERS_SRC_OPS_TPETRA_BLOCK_MULTI_VECTOR_PROD_VECTOR_HPP_

#include "../containers_ops_meta.hpp"
#include "../../multi_vector/containers_multi_vector_meta.hpp"

namespace pressio{ namespace containers{ namespace ops{

/*
 * multi_vector prod vector
 */


//-------------------------------------------------------//
//  TPETRA_BLOCK multivector with eigen vector
//-------------------------------------------------------//

// the result type is an Tpetra_Block wrapper and object is passed in
template <
  typename mvec_type,
  typename vec_type,
  typename res_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_tpetra_block<mvec_type>::value and
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value and
    containers::meta::is_vector_wrapper_eigen<vec_type>::value and
    containers::meta::is_vector_wrapper_tpetra_block<res_type>::value
    > * = nullptr
  >
void product(const mvec_type & mvA, const vec_type & vecB, res_type & C)
{
  /* computes: C = mvA*vecB */

  //zero out result
  C.setZero();

  // how many vectors are in mvA
  const auto numVecs = mvA.globalNumVectors();

  // size of vecB
  const size_t vecBLen = vecB.size();
  if (vecBLen != size_t(numVecs))
    assert(size_t(numVecs) == vecBLen);

  // get the wrapped trilinos tpetra multivector
  const auto mvA_mvv = mvA.data()->getMultiVectorView();
  const auto mvA_hv = mvA_mvv.template getLocalView<Kokkos::HostSpace>();

  // get data for the result too
  auto C_vv = C.data()->getVectorView();
  auto C_hv = C_vv.template getLocalView<Kokkos::HostSpace>();
  C_vv.template modify<Kokkos::HostSpace>();

  // my number of rows
  const auto myNrows = mvA_mvv.getLocalLength();
  // loop
  for (decltype(myNrows) i=0; i<myNrows; i++){
    for (decltype(numVecs) j=0; j<numVecs; j++){
      // we use C_hv(i,0) because C is Tpetra::Vector, which is
      // actually a Tpetra::MultiVector with one column, so C_hv
      // is a kokkos::View<scalar**,...> so we need to index the
      // zero column too
      C_hv(i,0) += mvA_hv(i,j) * vecB[j];
    }
  }
  using device_t = typename details::traits<res_type>::device_t;
  C.data()->template sync<device_t>();
}


// return a Tpetra vector wrapper, because a block one can be
// constructed from a regular tpetra vector
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
  Tpetra::Vector<typename details::traits<mvec_type>::scalar_t,
                 typename details::traits<mvec_type>::local_ordinal_t,
                 typename details::traits<mvec_type>::global_ordinal_t,
                 typename details::traits<mvec_type>::node_t>
                 >
{
  // here, mvA is distrubted, but vecB is NOT.
  // we interpret this as a linear combination of vectors

  // the data map of the multivector
  const auto rcpMap = mvA.getRCPDataMap();

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


}}}//end namespace pressio::containers::ops
#endif
#endif
