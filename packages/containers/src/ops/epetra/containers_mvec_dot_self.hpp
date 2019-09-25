
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_SRC_OPS_EPETRA_MULTI_VECTOR_DOT_SELF_HPP_
#define CONTAINERS_SRC_OPS_EPETRA_MULTI_VECTOR_DOT_SELF_HPP_

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
    ::pressio::containers::meta::is_multi_vector_wrapper_epetra<mvec_t>::value and
    ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<result_t>::value and
    ::pressio::containers::meta::wrapper_pair_have_same_scalar<mvec_t, result_t>::value
    > * = nullptr
  >
void dot_self(const mvec_t & A, result_t & C)
{
  // how many vectors are in A
  const auto numVecsA = A.globalNumVectors();
  const auto & Adata = *A.data();
  assert(C.rows() == numVecsA);
  assert(C.cols() == numVecsA);

  // A dot A = A^T*A, which yields a symmetric matrix
  // only need to compute half and fill remaining entries accordingly
  for (auto i=0; i<numVecsA; i++){
    for (auto j=i; j<numVecsA; j++){
      Adata(i)->Dot( *(Adata(j)), &C(i,j) );
      // fill the lower triangular part
      C(j,i) = C(i,j);
    }
  }
}

// return result in the form of eigen dense matrix wrapper
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
