
#ifndef CORE_MATRIX_OPERATIONS_SPARSE_DENSE_MATRIX_PRODUCT_HPP_
#define CORE_MATRIX_OPERATIONS_SPARSE_DENSE_MATRIX_PRODUCT_HPP_

#include "../core_ops_meta.hpp"
#include "../../matrix/core_matrix_meta.hpp"

namespace rompp{ namespace core{ namespace ops{
  
/*---------------------------------------------------------
C = A B
- A is sparse matrix from eigen 
- B is dense matrix from eigen 
- result C can be either dense or sparse
---------------------------------------------------------*/

// C is a DENSE or SPARSE matrix 
template <typename A_t, typename B_t, typename C_t,
  core::meta::enable_if_t<
    core::meta::is_eigen_sparse_matrix_wrapper<A_t>::value and 
    core::meta::is_eigen_dense_matrix_wrapper<B_t>::value and 
    core::meta::wrapper_triplet_have_same_scalar<A_t,B_t,C_t>::value and 
    (core::meta::is_eigen_dense_matrix_wrapper<C_t>::value or 
     core::meta::is_eigen_sparse_matrix_wrapper<C_t>::value )
    > * = nullptr
  >
void product(const A_t & A,
	     const B_t & B,
	     C_t & C){

  assert(C.rows() == A.rows());
  assert(C.cols() == B.cols());
  (*C.data()) = (*A.data()) * (*B.data());
}

      
template <typename A_t, typename B_t, 
    core::meta::enable_if_t<
      core::meta::is_eigen_sparse_matrix_wrapper<A_t>::value and 
      core::meta::is_eigen_dense_matrix_wrapper<B_t>::value and
      core::meta::wrapper_pair_have_same_scalar<A_t,B_t>::value
      > * = nullptr
    >
B_t product(const A_t & A,
	     const B_t & B){
  B_t C(A.rows(),B.cols());
  product(A, B, C);
  return C;
}



}}} // end namespace rompp::core::ops
#endif
