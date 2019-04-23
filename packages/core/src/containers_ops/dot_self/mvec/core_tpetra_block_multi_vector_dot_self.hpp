
#ifdef HAVE_TRILINOS
#ifndef CORE_TPETRA_BLOCK_MULTI_VECTOR_DOT_SELF_HPP_
#define CORE_TPETRA_BLOCK_MULTI_VECTOR_DOT_SELF_HPP_

#include "../../core_ops_meta.hpp"
#include "../../../multi_vector/core_multi_vector_meta.hpp"
#include "core_tpetra_multi_vector_dot_self.hpp"

namespace rompp{ namespace core{ namespace ops{

// Tpetra multivector (A) dot self
// this is equivalent to doing A^T * A

template <
  typename mvec_t,
  ::rompp::mpl::enable_if_t<
    core::meta::is_multi_vector_wrapper_tpetra_block<
      mvec_t
      >::value
    > * = nullptr
  >
void dot_self(const mvec_t & mvA,
	      core::Matrix<
	       Eigen::Matrix<
		typename core::details::traits<mvec_t>::scalar_t,
		Eigen::Dynamic,
		Eigen::Dynamic
	       >
	      > & C)
{
  // get a tpetra multivector that views the data
  const auto mvView = mvA.getMultiVectorView();

  // use the dot_self already implemented for tpetra MV
  ::rompp::core::ops::dot_self(mvView, C);
}


template <
  typename mvec_t,
  ::rompp::mpl::enable_if_t<
    core::meta::is_multi_vector_wrapper_tpetra_block<
      mvec_t
      >::value
    > * = nullptr
  >
auto dot_self(const mvec_t & mvA)
  -> core::Matrix<
  Eigen::Matrix<typename core::details::traits<mvec_t>::scalar_t,
  Eigen::Dynamic, Eigen::Dynamic>>{

  using sc_t	= typename core::details::traits<mvec_t>::scalar_t;
  using eig_mat = Eigen::Matrix< sc_t, Eigen::Dynamic, Eigen::Dynamic>;
  using res_t	= core::Matrix<eig_mat>;

  auto numVecsA = mvA.globalNumVectors();
  res_t C(numVecsA, numVecsA);
  dot_self(mvA, C);
  return C;
}


}}} // end namespace rompp::core::ops
#endif
#endif
