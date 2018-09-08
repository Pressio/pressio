
#ifndef CORE_MATRIX_OPERATIONS_SPARSE_MATRIX_VECTOR_PRODUCT_HPP_
#define CORE_MATRIX_OPERATIONS_SPARSE_MATRIX_VECTOR_PRODUCT_HPP_

#include "../../meta/core_vector_meta.hpp"
#include "../../vector/concrete/core_vector_sharedmem_eigen.hpp"
#include "../../meta/core_matrix_meta.hpp"
#include "../concrete/core_matrix_dense_sharedmem_eigen.hpp"
#include "../concrete/core_matrix_sparse_sharedmem_eigen.hpp"
#include "Epetra_LocalMap.h"
#include <Epetra_Export.h>
#include <Epetra_Import.h>

namespace core{
namespace mat_ops{
  
/*---------------------------------------------------------
  EPETRA 
  c = A b , 
  - A = crs matrix 
  - b = SINGLE vector
  - c is epetra vector
-----------------------------------------------------------*/

template <typename matrix_type,
	  typename vector_type,
	  core::meta::enable_if_t<
	    details::traits<matrix_type>::isMatrix==1 &&
	    details::traits<matrix_type>::isEpetra==1 &&
	    details::traits<matrix_type>::isSparse==1 &&
	    details::traits<vector_type>::isVector==1 &&
	    details::traits<vector_type>::isEpetra==1
	    > * = nullptr
	  >
auto product(const matrix_type & A,
	     const vector_type & b,
	     bool transposeA = false){

  assert( A.isFillingCompleted() );
  assert( A.globalCols() == b.globalSize() );
  vector_type c( A.getRangeDataMap() );
  A.data()->Multiply(transposeA, *b.data(), *c.data());
  return c;
}


template <typename matrix_type,
	  typename vector_type,
	  core::meta::enable_if_t<
	    details::traits<matrix_type>::isMatrix==1 &&
	    details::traits<matrix_type>::isEpetra==1 &&
	    details::traits<matrix_type>::isSparse==1 &&
	    details::traits<vector_type>::isVector==1 &&
	    details::traits<vector_type>::isEpetra==1
	    > * = nullptr
	  >
void product(const matrix_type & A,
	     const vector_type & b,
	     vector_type & c,
	     bool transposeA = false){

  assert( A.isFillingCompleted() );
  assert( A.globalCols() == b.globalSize() );
  A.data()->Multiply(transposeA, *b.data(), *c.data());
}

  
/*---------------------------------------------------------
c = A b
- A is matrix from eigen
- b is vector from eigen
- c is an eigen vector storing the result 
---------------------------------------------------------*/
  
template <typename matrix_type,
	  typename vector_t,
	  typename res_t,
	  core::meta::enable_if_t<
	    details::traits<vector_t>::isEigen &&
	    details::traits<vector_t>::isVector &&
	    details::traits<res_t>::isEigen &&
	    details::traits<res_t>::isVector &&
	    details::traits<matrix_type>::isEigen &&
	    details::traits<matrix_type>::isSparse &&
	    std::is_same<typename details::traits<matrix_type>::scalar_t,
			 typename details::traits<vector_t>::scalar_t
			 >::value && 
	    std::is_same<typename details::traits<matrix_type>::scalar_t,
			 typename details::traits<res_t>::scalar_t
			 >::value 
	    > * = nullptr
	  >
void product(const matrix_type & A,
	     const vector_t & b,
	     res_t & c){

  assert(A.cols() == b.size());
  assert(c.size() == A.rows());
  (*c.data()) = (*A.data()) * (*b.data());
}


  
} // end namespace mat_ops
} // end namespace core
#endif
