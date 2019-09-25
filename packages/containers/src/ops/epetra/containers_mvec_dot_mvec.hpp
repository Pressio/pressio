
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_SRC_OPS_EPETRA_MULTI_VECTOR_DOT_MULTI_VECTOR_HPP_
#define CONTAINERS_SRC_OPS_EPETRA_MULTI_VECTOR_DOT_MULTI_VECTOR_HPP_

#include "../containers_ops_meta.hpp"
#include "../../multi_vector/containers_multi_vector_meta.hpp"

namespace pressio{ namespace containers{ namespace ops{

/*
 * multi_vector dot multi vector
 */


// result is stored into an Eigen dense matrix wrapper
template <
  typename mvec_t,
  typename result_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_epetra<mvec_t>::value and
    ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<result_t>::value and
    ::pressio::containers::meta::wrapper_pair_have_same_scalar<mvec_t, result_t>::value
    > * = nullptr
  >
void dot(const mvec_t & mvA, const mvec_t & mvB, result_t & C){

  // how many vectors are in mvA and mvB
  const auto numVecsA = mvA.globalNumVectors();
  const auto numVecsB = mvB.globalNumVectors();
  assert( mvA.globalLength() == mvB.globalLength());
  auto const & mvAdata = *mvA.data();
  auto const & mvBdata = *mvB.data();

  assert(C.rows() == numVecsA);
  assert(C.cols() == numVecsB);
  // compute dot between every column of A with every col of B
  for (size_t i=0; i<(size_t)numVecsA; i++){
    for (size_t j=0; j<(size_t)numVecsB; j++){
      mvAdata(i)->Dot( *(mvBdata(j)), &C(i,j) );
    }
  }
}

// result is stored into an Eigen dense DYNAMIC matrix wrapper
template <
  typename mvec_t,
  typename result_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_epetra<mvec_t>::value and
    ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<result_t>::value and
    ::pressio::containers::details::traits<result_t>::is_dynamic and
    ::pressio::containers::meta::wrapper_pair_have_same_scalar<mvec_t, result_t>::value
    > * = nullptr
  >
result_t dot(const mvec_t & mvA, const mvec_t & mvB){
  const auto numVecsA = mvA.globalNumVectors();
  const auto numVecsB = mvB.globalNumVectors();
  result_t C(numVecsA, numVecsB);
  dot(mvA, mvB, C);
  return C;
}

}}}//end namespace pressio::containers::ops
#endif
#endif
