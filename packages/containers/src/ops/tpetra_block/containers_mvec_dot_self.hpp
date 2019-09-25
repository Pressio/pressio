
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_SRC_OPS_TPETRA_BLOCK_MULTI_VECTOR_DOT_SELF_HPP_
#define CONTAINERS_SRC_OPS_TPETRA_BLOCK_MULTI_VECTOR_DOT_SELF_HPP_

#include "../containers_ops_meta.hpp"
#include "../../multi_vector/containers_multi_vector_meta.hpp"

namespace pressio{ namespace containers{ namespace ops{

/*
 * A dot A
 * multi_vector dot self: equivalent to doing A^T A
 */

// store result into eigen dense matrix wrapper
template <
  typename mvec_t,
  typename result_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_tpetra_block<mvec_t>::value and
    ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<result_t>::value and
    ::pressio::containers::meta::wrapper_pair_have_same_scalar<mvec_t, result_t>::value
    > * = nullptr
  >
void dot_self(const mvec_t & mvA, result_t & C)
{
  // get a tpetra multivector that views the data
  const auto mvView = mvA.data()->getMultiVectorView();

  // how many vectors are in mvA and mvB
  const auto numVecsA = mvA.globalNumVectors();

  // A dot A = A^T*A, which yields a symmetric matrix
  // only need to compute half and fill remaining entries accordingly
  for (auto i=0; i<numVecsA; i++)
  {
    // colI is a Teuchos::RCP<Vector<...>>
    const auto colI = mvView.getVector(i);
    for (auto j=i; j<numVecsA; j++)
    {
      const auto colJ = mvView.getVector(j);
      C(i,j) = colI->dot(*colJ);
      C(j,i) = C(i,j);
    }
  }
}

// return result in the form of eigen dense matrix wrapper
template <
  typename mvec_t,
  typename result_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_tpetra_block<mvec_t>::value and
    ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<result_t>::value and
    ::pressio::containers::details::traits<result_t>::is_dynamic and
    ::pressio::containers::meta::wrapper_pair_have_same_scalar<mvec_t, result_t>::value
    > * = nullptr
  >
result_t dot_self(const mvec_t & mvA)
{
  const auto numVecsA = mvA.globalNumVectors();
  result_t C(numVecsA, numVecsA);
  dot_self(mvA, C);
  return C;
}


}}}//end namespace pressio::containers::ops
#endif
#endif
