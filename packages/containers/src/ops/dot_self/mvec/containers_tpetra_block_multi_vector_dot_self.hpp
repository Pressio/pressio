
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_TPETRA_BLOCK_MULTI_VECTOR_DOT_SELF_HPP_
#define CONTAINERS_TPETRA_BLOCK_MULTI_VECTOR_DOT_SELF_HPP_

#include "../../containers_ops_meta.hpp"
#include "../../../multi_vector/containers_multi_vector_meta.hpp"
#include "containers_tpetra_multi_vector_dot_self.hpp"

namespace pressio{ namespace containers{ namespace ops{

// Tpetra multivector (A) dot self
// this is equivalent to doing A^T * A

template <
  typename mvec_t,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_tpetra_block<
      mvec_t
      >::value
    > * = nullptr
  >
void dot_self(const mvec_t & mvA,
	      containers::Matrix<
	       Eigen::Matrix<
		typename containers::details::traits<mvec_t>::scalar_t,
		Eigen::Dynamic,
		Eigen::Dynamic
	       >
	      > & C)
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


template <
  typename mvec_t,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_tpetra_block<
      mvec_t
      >::value
    > * = nullptr
  >
auto dot_self(const mvec_t & mvA)
  -> containers::Matrix<
  Eigen::Matrix<typename containers::details::traits<mvec_t>::scalar_t,
  Eigen::Dynamic, Eigen::Dynamic>>{

  using sc_t	= typename containers::details::traits<mvec_t>::scalar_t;
  using eig_mat = Eigen::Matrix< sc_t, Eigen::Dynamic, Eigen::Dynamic>;
  using res_t	= containers::Matrix<eig_mat>;

  auto numVecsA = mvA.globalNumVectors();
  res_t C(numVecsA, numVecsA);
  dot_self(mvA, C);
  return C;
}


}}} // end namespace pressio::containers::ops
#endif
#endif
