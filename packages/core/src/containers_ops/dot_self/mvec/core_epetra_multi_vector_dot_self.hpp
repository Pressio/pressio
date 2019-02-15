
#ifdef HAVE_TRILINOS
#ifndef CORE_EPETRA_MULTI_VECTOR_DOT_SELF_HPP_
#define CORE_EPETRA_MULTI_VECTOR_DOT_SELF_HPP_

#include "../../core_ops_meta.hpp"
#include "../../../multi_vector/core_multi_vector_meta.hpp"

namespace rompp{ namespace core{ namespace ops{

// Epetra multivector (A) dot self
// this is equivalent to doing A^T * A

template <typename mvec_t,
  core::meta::enable_if_t<
    core::meta::is_multi_vector_wrapper_epetra<mvec_t>::value
    > * = nullptr
  >
void dot_self(const mvec_t & A,
	      core::Matrix<
	      Eigen::Matrix<typename core::details::traits<mvec_t>::scalar_t,
	      Eigen::Dynamic, Eigen::Dynamic>
	      > & C)
{
  // how many vectors are in A
  auto numVecsA = A.globalNumVectors();
  auto const & Adata = *A.data();
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


template <typename mvec_t,
  core::meta::enable_if_t<
    core::meta::is_multi_vector_wrapper_epetra<mvec_t>::value
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
