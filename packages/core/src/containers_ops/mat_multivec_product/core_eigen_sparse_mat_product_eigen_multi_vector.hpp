
#ifndef CORE_CONTAINERS_OPS_EIGEN_SPARSE_MAT_PRODUCT_EIGEN_MULTI_VECTOR_HPP_
#define CORE_CONTAINERS_OPS_EIGEN_SPARSE_MAT_PRODUCT_EIGEN_MULTI_VECTOR_HPP_

#include "../core_ops_meta.hpp"
#include "../../vector/core_vector_meta.hpp"
#include "../../multi_vector/core_multi_vector_meta.hpp"
#include "../../matrix/core_matrix_meta.hpp"
#include "../../matrix/concrete/core_matrix_sparse_sharedmem_eigen.hpp"
#include "../../matrix/concrete/core_matrix_dense_sharedmem_eigen_dynamic.hpp"
#include "../../multi_vector/concrete/core_multi_vector_sharedmem_eigen_dynamic.hpp"

namespace rompp{ namespace core{ namespace ops{

/*---------------------------------------------
 * C = A * B
 * A : eigen sparse matrix wrapper
 * B: eigen multivector wrapper
 * C: eigen dense matrix wrapper
 *-----------------------------------------------*/
template <typename mat_type,
	  typename mvec_type,
  ::rompp::mpl::enable_if_t<
   core::meta::is_sparse_matrix_wrapper_eigen<mat_type>::value and
   core::meta::is_multi_vector_wrapper_eigen<mvec_type>::value and
    core::meta::wrapper_pair_have_same_scalar<mat_type, mvec_type>::value
    > * = nullptr
  >
void product(const mat_type & A,
	     const mvec_type & mv,
	     core::Matrix<Eigen::MatrixXd> & C){

  assert( C.rows() == A.rows() );
  assert( mv.length() == A.cols() );
  assert( C.cols() == mv.numVectors() );
  (*C.data()) = (*A.data()) * (*mv.data());
}//end function

// construct and return result
template <typename mat_type,
	  typename mvec_type,
  ::rompp::mpl::enable_if_t<
    core::meta::is_sparse_matrix_wrapper_eigen<mat_type>::value and
    core::meta::is_multi_vector_wrapper_eigen<mvec_type>::value and
    core::meta::wrapper_pair_have_same_scalar<mat_type, mvec_type>::value
    > * = nullptr
  >
auto product(const mat_type & A, const mvec_type & mv)
  -> core::Matrix<Eigen::MatrixXd>
{
  core::Matrix<Eigen::MatrixXd> C(A.rows(), mv.numVectors());
  product(A,mv,C);
  return C;
}//end function


}}}//end namespace rompp::core::ops
#endif
