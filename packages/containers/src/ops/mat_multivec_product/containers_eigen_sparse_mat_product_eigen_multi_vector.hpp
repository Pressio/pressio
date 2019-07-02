
#ifndef CONTAINERS_CONTAINERS_OPS_EIGEN_SPARSE_MAT_PRODUCT_EIGEN_MULTI_VECTOR_HPP_
#define CONTAINERS_CONTAINERS_OPS_EIGEN_SPARSE_MAT_PRODUCT_EIGEN_MULTI_VECTOR_HPP_

#include "../containers_ops_meta.hpp"
#include "../../vector/containers_vector_meta.hpp"
#include "../../multi_vector/containers_multi_vector_meta.hpp"
#include "../../matrix/containers_matrix_meta.hpp"
#include "../../matrix/concrete/containers_matrix_sparse_sharedmem_eigen.hpp"
#include "../../matrix/concrete/containers_matrix_dense_sharedmem_eigen_dynamic.hpp"
#include "../../multi_vector/concrete/containers_multi_vector_sharedmem_eigen_dynamic.hpp"

namespace pressio{ namespace containers{ namespace ops{

/*---------------------------------------------
 * C = A * B
 * A : eigen sparse matrix wrapper
 * B: eigen multivector wrapper
 * C: eigen dense matrix wrapper
 *-----------------------------------------------*/
template <typename mat_type,
	  typename mvec_type,
  ::pressio::mpl::enable_if_t<
   ::pressio::containers::meta::is_sparse_matrix_wrapper_eigen<mat_type>::value and
   ::pressio::containers::meta::is_multi_vector_wrapper_eigen<mvec_type>::value and
    ::pressio::containers::meta::wrapper_pair_have_same_scalar<mat_type, mvec_type>::value
    > * = nullptr
  >
void product(const mat_type & A,
	     const mvec_type & mv,
	     ::pressio::containers::Matrix<Eigen::MatrixXd> & C){

  assert( C.rows() == A.rows() );
  assert( mv.length() == A.cols() );
  assert( C.cols() == mv.numVectors() );
  (*C.data()) = (*A.data()) * (*mv.data());
}//end function

// construct and return result
template <typename mat_type,
	  typename mvec_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_sparse_matrix_wrapper_eigen<mat_type>::value and
    ::pressio::containers::meta::is_multi_vector_wrapper_eigen<mvec_type>::value and
    ::pressio::containers::meta::wrapper_pair_have_same_scalar<mat_type, mvec_type>::value
    > * = nullptr
  >
auto product(const mat_type & A, const mvec_type & mv)
  -> ::pressio::containers::Matrix<Eigen::MatrixXd>
{
  ::pressio::containers::Matrix<Eigen::MatrixXd> C(A.rows(), mv.numVectors());
  product(A,mv,C);
  return C;
}//end function


}}}//end namespace pressio::containers::ops
#endif
