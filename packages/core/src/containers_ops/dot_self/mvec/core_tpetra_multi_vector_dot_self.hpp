
#ifdef HAVE_TRILINOS
#ifndef CORE_TPETRA_MULTI_VECTOR_DOT_SELF_HPP_
#define CORE_TPETRA_MULTI_VECTOR_DOT_SELF_HPP_

#include "../../core_ops_meta.hpp"
#include "../../../multi_vector/core_multi_vector_meta.hpp"

namespace rompp{ namespace core{ namespace ops{

// Tpetra multivector (A) dot self
// this is equivalent to doing A^T * A

template <typename mvec_t,
  ::rompp::mpl::enable_if_t<
    core::meta::is_multi_vector_wrapper_tpetra<mvec_t>::value
    > * = nullptr
  >
void dot_self(const mvec_t & mvA,
	      core::Matrix<
	      Eigen::Matrix<typename core::details::traits<mvec_t>::scalar_t,
	      Eigen::Dynamic, Eigen::Dynamic>
	      > & C)
{
  // how many vectors are in mvA and mvB
  auto numVecsA = mvA.globalNumVectors();
  assert(C.rows() == numVecsA);
  assert(C.cols() == numVecsA);

  // A dot A = A^T*A, which yields a symmetric matrix
  // only need to compute half and fill remaining entries accordingly
  for (auto i=0; i<numVecsA; i++){
    // colI is a Teuchos::RCP<Vector<...>>
    auto colI = mvA.data()->getVector(i);
    for (auto j=i; j<numVecsA; j++){
      auto colJ = mvA.data()->getVector(j);
      C(i,j) = colI->dot(*colJ);
      C(j,i) = C(i,j);
    }
  }
}


template <typename mvec_t,
  ::rompp::mpl::enable_if_t<
    core::meta::is_multi_vector_wrapper_tpetra<mvec_t>::value
    > * = nullptr
  >
auto dot_self(const mvec_t & mvA)
  -> core::Matrix<
  Eigen::Matrix<typename core::details::traits<mvec_t>::scalar_t,
  Eigen::Dynamic, Eigen::Dynamic>>{

  using sc_t = typename core::details::traits<mvec_t>::scalar_t;
  using eig_mat = Eigen::Matrix< sc_t, Eigen::Dynamic, Eigen::Dynamic>;
  using res_t = core::Matrix<eig_mat>;

  auto numVecsA = mvA.globalNumVectors();
  res_t C(numVecsA, numVecsA);
  dot_self(mvA, C);
  return C;
}


}}} // end namespace rompp::core::ops
#endif
#endif
