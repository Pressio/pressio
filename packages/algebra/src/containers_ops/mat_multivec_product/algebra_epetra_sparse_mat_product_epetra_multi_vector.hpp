
#ifdef HAVE_TRILINOS
#ifndef ALGEBRA_CONTAINERS_OPS_EPETRA_SPARSE_MAT_PRODUCT_EPETRA_MULTI_VECTOR_HPP_
#define ALGEBRA_CONTAINERS_OPS_EPETRA_SPARSE_MAT_PRODUCT_EPETRA_MULTI_VECTOR_HPP_

#include "../algebra_ops_meta.hpp"
#include "../../vector/algebra_vector_meta.hpp"
#include "../../multi_vector/algebra_multi_vector_meta.hpp"
#include "../../matrix/algebra_matrix_meta.hpp"
#include "../../matrix/concrete/algebra_matrix_sparse_sharedmem_eigen.hpp"
#include "../../matrix/concrete/algebra_matrix_dense_sharedmem_eigen_dynamic.hpp"
#include "../../multi_vector/concrete/algebra_multi_vector_sharedmem_eigen_dynamic.hpp"


namespace rompp{ namespace algebra{ namespace ops{

/*------------------------------------------
 *  C = A B
 *
 * A: EPETRA CRS matrix wrapper
 * B: epetra multivector wrapper
 * C: epetra multivector wrapper
 *------------------------------------------*/

template <typename mat_type,
	  typename mvec_type,
  ::rompp::mpl::enable_if_t<
    ::rompp::algebra::meta::is_sparse_matrix_wrapper_epetra<mat_type>::value and
    ::rompp::algebra::meta::is_multi_vector_wrapper_epetra<mvec_type>::value and
    ::rompp::algebra::meta::wrapper_pair_have_same_scalar<mat_type, mvec_type>::value
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
  ::rompp::mpl::enable_if_t<
    ::rompp::algebra::meta::is_sparse_matrix_wrapper_epetra<mat_type>::value and
    ::rompp::algebra::meta::is_multi_vector_wrapper_epetra<mvec_type>::value and
    ::rompp::algebra::meta::wrapper_pair_have_same_scalar<mat_type, mvec_type>::value
    > * = nullptr
  >
mvec_type product(const mat_type & A, const mvec_type & B)
{

  mvec_type C( A.getRangeDataMap(), B.globalNumVectors() );
  product(A,B,C);
  return C;
}


}}}//end namespace rompp::ops
#endif
#endif

