
#ifndef ALGEBRA_EIGEN_MULTI_VECTOR_DOT_SELF_HPP_
#define ALGEBRA_EIGEN_MULTI_VECTOR_DOT_SELF_HPP_

#include "../../algebra_ops_meta.hpp"
#include "../../../multi_vector/algebra_multi_vector_meta.hpp"

namespace rompp{ namespace algebra{ namespace ops{

// Eigen multivector (A) dot self
// this is equivalent to doing A^T * A

template <
  typename mvec_t,
  ::rompp::mpl::enable_if_t<
    algebra::meta::is_multi_vector_wrapper_eigen<mvec_t>::value
    > * = nullptr
  >
void dot_self(const mvec_t & A,
	      algebra::Matrix<
	      Eigen::Matrix<typename algebra::details::traits<mvec_t>::scalar_t,
	      Eigen::Dynamic, Eigen::Dynamic>
	      > & C)
{
  // // how many vectors are in A
  // auto numVecsA = A.numVectors();
  // // auto const & Adata = *A.data();
  // assert(C.rows() == numVecsA);
  // assert(C.cols() == numVecsA);

  *C.data() = A.data()->transpose() * (*A.data());
}


template <
  typename mvec_t,
  ::rompp::mpl::enable_if_t<
    algebra::meta::is_multi_vector_wrapper_eigen<mvec_t>::value
    > * = nullptr
  >
auto dot_self(const mvec_t & mvA)
  -> algebra::Matrix<
    Eigen::Matrix<typename algebra::details::traits<mvec_t>::scalar_t,
		  Eigen::Dynamic, Eigen::Dynamic>
    >{

  using sc_t = typename algebra::details::traits<mvec_t>::scalar_t;
  using eig_mat = Eigen::Matrix< sc_t, Eigen::Dynamic, Eigen::Dynamic>;
  using res_t = algebra::Matrix<eig_mat>;

  auto numVecsA = mvA.numVectors();
  res_t C(numVecsA, numVecsA);
  dot_self(mvA, C);
  return C;
}

}}} // end namespace rompp::algebra::ops
#endif
