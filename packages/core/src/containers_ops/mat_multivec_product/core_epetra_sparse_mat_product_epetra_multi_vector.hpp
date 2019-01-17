
#ifdef HAVE_TRILINOS
#ifndef CORE_CONTAINERS_OPS_EPETRA_SPARSE_MAT_PRODUCT_EPETRA_MULTI_VECTOR_HPP_
#define CORE_CONTAINERS_OPS_EPETRA_SPARSE_MAT_PRODUCT_EPETRA_MULTI_VECTOR_HPP_

#include "../core_ops_meta.hpp"
#include "../../vector/core_vector_meta.hpp"
#include "../../multi_vector/core_multi_vector_meta.hpp"
#include "../../matrix/core_matrix_meta.hpp"
#include "../../matrix/concrete/core_matrix_sparse_sharedmem_eigen.hpp"
#include "../../matrix/concrete/core_matrix_dense_sharedmem_eigen_dynamic.hpp"
#include "../../multi_vector/concrete/core_multi_vector_sharedmem_eigen_dynamic.hpp"


namespace rompp{ namespace core{ namespace ops{

/*------------------------------------------
 *  C = A B
 *
 * A: EPETRA CRS matrix wrapper
 * B: epetra multivector wrapper
 * C: epetra multivector wrapper
 *------------------------------------------*/

template <typename mat_type,
	  typename mvec_type,
  core::meta::enable_if_t<
    core::meta::is_sparse_matrix_wrapper_epetra<mat_type>::value and
    core::meta::is_multi_vector_wrapper_epetra<mvec_type>::value and
    core::meta::wrapper_pair_have_same_scalar<mat_type, mvec_type>::value
    > * = nullptr
  >
void product(const mat_type & A,
	     const mvec_type & B,
	     mvec_type & C){

  assert( A.isFillingCompleted() );
  assert( A.globalCols() == B.globalLength() );
  assert( C.globalNumVectors() == B.globalNumVectors() );
  assert( C.globalLength() == A.getRangeDataMap().NumGlobalElements() );
  A.data()->Multiply(false, *B.data(), *C.data());
}

template <typename mat_type,
	  typename mvec_type,
  core::meta::enable_if_t<
    core::meta::is_sparse_matrix_wrapper_epetra<mat_type>::value and
    core::meta::is_multi_vector_wrapper_epetra<mvec_type>::value and
    core::meta::wrapper_pair_have_same_scalar<mat_type, mvec_type>::value
    > * = nullptr
  >
mvec_type product(const mat_type & A, const mvec_type & B)
{

  mvec_type C( A.getRangeDataMap(), B.globalNumVectors() );
  product(A,B,C);
  return C;
}


}}}//end namespace rompp::core::ops
#endif
#endif

