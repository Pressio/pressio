
#ifndef CORE_EPETRA_MULTI_VECTOR_DOT_EPETRA_MULTI_VECTOR_HPP_
#define CORE_EPETRA_MULTI_VECTOR_DOT_EPETRA_MULTI_VECTOR_HPP_

#include "../../core_ops_meta.hpp"
#include "../../../multi_vector/core_multi_vector_meta.hpp"

namespace rompp{ namespace core{ namespace ops{

// Epetra multivector dot epetra multi vector

#ifdef HAVE_TRILINOS
template <typename mvec_t,
  core::meta::enable_if_t<
    core::meta::is_epetra_multi_vector_wrapper<mvec_t>::value
    > * = nullptr
  >
auto dot(const mvec_t & mvA, const mvec_t & mvB)
  -> core::Matrix<
  Eigen::Matrix<typename core::details::traits<mvec_t>::scalar_t,
  Eigen::Dynamic, Eigen::Dynamic>>{

  using sc_t = typename core::details::traits<mvec_t>::scalar_t;
  using eig_mat = Eigen::Matrix< sc_t, Eigen::Dynamic, Eigen::Dynamic>;
  using res_t = core::Matrix<eig_mat>;

  // how many vectors are in mvA and mvB
  auto numVecsA = mvA.globalNumVectors();
  auto nRowsA = mvA.globalLength();
  auto numVecsB = mvB.globalNumVectors();
  auto nRowsB = mvB.globalLength();
  assert( nRowsA == nRowsB );
  auto const & mvAdata = *mvA.data();
  auto const & mvBdata = *mvB.data();

  // result
  res_t C(numVecsA, numVecsB);
  // compute dot between every column of A with every col of B
  for (auto i=0; i<numVecsA; i++){
    for (auto j=0; j<numVecsB; j++){
      mvAdata(i)->Dot( *(mvBdata(j)), &C(i,j) );
    }
  }
  return C;
}
#endif
//--------------------------------------------------------
      
  
}}} // end namespace rompp::core::ops
#endif
