
#ifdef HAVE_TRILINOS
#ifndef ALGEBRA_TPETRA_BLOCK_MULTI_VECTOR_DOT_TPETRA_BLOCK_MULTI_VECTOR_HPP_
#define ALGEBRA_TPETRA_BLOCK_MULTI_VECTOR_DOT_TPETRA_BLOCK_MULTI_VECTOR_HPP_

#include "../../algebra_ops_meta.hpp"
#include "../../../multi_vector/algebra_multi_vector_meta.hpp"

namespace rompp{ namespace algebra{ namespace ops{

template <typename mvec_t,
  ::rompp::mpl::enable_if_t<
    algebra::meta::is_multi_vector_wrapper_tpetra_block<mvec_t>::value
    > * = nullptr
  >
void dot(const mvec_t & mvA,
	 const mvec_t & mvB,
	 algebra::Matrix<
	 Eigen::Matrix<typename algebra::details::traits<mvec_t>::scalar_t,
	 Eigen::Dynamic, Eigen::Dynamic>
	 > & C)
{
  // how many vectors are in mvA and mvB
  auto numVecsA = mvA.globalNumVectors();
  auto numVecsB = mvB.globalNumVectors();

  auto mvA_v = mvA.data()->getMultiVectorView();
  auto mvB_v = mvB.data()->getMultiVectorView();

  // compute dot between every column of A with every col of B
  for (auto i=0; i<numVecsA; i++)
  {
    // colI is a Teuchos::RCP<Vector<...>>
    auto colI = mvA_v.getVector(i);
    for (auto j=0; j<numVecsB; j++)
    {
      auto colJ = mvB_v.getVector(j);
      C(i,j) = colI->dot(*colJ);
    }
  }
}


template <typename mvec_t,
  ::rompp::mpl::enable_if_t<
    algebra::meta::is_multi_vector_wrapper_tpetra_block<mvec_t>::value
    > * = nullptr
  >
auto dot(const mvec_t & mvA, const mvec_t & mvB)
  -> algebra::Matrix<
  Eigen::Matrix<typename algebra::details::traits<mvec_t>::scalar_t,
  Eigen::Dynamic, Eigen::Dynamic>>{

  using sc_t = typename algebra::details::traits<mvec_t>::scalar_t;
  using eig_mat = Eigen::Matrix< sc_t, Eigen::Dynamic, Eigen::Dynamic>;
  using res_t = algebra::Matrix<eig_mat>;

  auto numVecsA = mvA.globalNumVectors();
  auto numVecsB = mvB.globalNumVectors();
  res_t C(numVecsA, numVecsB);
  dot(mvA, mvB, C);

  return C;
}
//--------------------------------------------------------


}}} // end namespace rompp::algebra::ops
#endif
#endif
