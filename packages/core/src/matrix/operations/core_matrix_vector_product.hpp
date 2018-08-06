
#ifndef CORE_MATRIX_OPERATIONS_MATRIX_VECTOR_PRODUCT_HPP_
#define CORE_MATRIX_OPERATIONS_MATRIX_VECTOR_PRODUCT_HPP_

#include "../../meta/core_vector_meta.hpp"
#include "../../vector/concrete/core_vector_serial_eigen.hpp"
#include "../../meta/core_matrix_meta.hpp"
#include "../concrete/core_matrix_dense_serial_eigen.hpp"
#include "../concrete/core_matrix_sparse_serial_eigen.hpp"

namespace core{

/*-----------------------------------------------------
c = A b
- A is matrix from eigen
- b is vector from eigen
-----------------------------------------------------*/
  
template <typename matrix_type,
	  typename vector_t_1,
	  typename vector_t_2,
	  typename std::enable_if<
	    // b must be a vector and from eigen
	    details::traits<vector_t_1>::isEigen &&
	    details::traits<vector_t_1>::isVector &&
	    // c must be a vector and from eigen
	    details::traits<vector_t_2>::isEigen &&
	    details::traits<vector_t_2>::isVector &&
	    // A must be a matrix and from eigen
	    details::traits<matrix_type>::isEigen &&
	    (details::traits<matrix_type>::isSparse ||
	     details::traits<matrix_type>::isDense ) &&
	    // we need to have matching scalar types
	    std::is_same<typename details::traits<matrix_type>::scalar_t,
			 typename details::traits<vector_t_1>::scalar_t
			 >::value &&
	    std::is_same<typename details::traits<matrix_type>::scalar_t,
			 typename details::traits<vector_t_2>::scalar_t
			 >::value 			 
	    >::type * = nullptr
	  >
void matrixVectorProduct(const matrix_type & A,
			 const vector_t_1 & b,
			 vector_t_2 & c)
{

  assert(A.cols() == b.size());
  assert(c.size() == A.rows());
  (*c.data()) = (*A.data()) * (*b.data());
}




  
/*-----------------------------------------------------------
  EPETRA 
  c = A b , 
  - A = crs matrix 
  - b = SINGLE vector
-----------------------------------------------------------*/

template <typename matrix_type,
	  typename vector_type,
	  typename std::enable_if<
	    details::traits<matrix_type>::isMatrix==1 &&
	    details::traits<matrix_type>::isEpetra==1 &&
	    details::traits<matrix_type>::isSparse==1 &&
	    details::traits<vector_type>::isVector==1 &&
	    details::traits<vector_type>::isEpetra==1
	    >::type * = nullptr
	  >
auto matrixVectorProduct(const matrix_type & A,
			 const vector_type & b,
			 bool transposeA = false)
{
  assert( A.isFillingCompleted() );
  assert( A.globalCols() == b.globalSize() );
  vector_type c( A.getRangeDataMap() );
  A.data()->Multiply(transposeA, *b.data(), *c.data());
  return c;
}


template <typename matrix_type,
	  typename vector_type,
	  typename std::enable_if<
	    details::traits<matrix_type>::isMatrix==1 &&
	    details::traits<matrix_type>::isEpetra==1 &&
	    details::traits<matrix_type>::isSparse==1 &&
	    details::traits<vector_type>::isVector==1 &&
	    details::traits<vector_type>::isEpetra==1
	    >::type * = nullptr
	  >
void matrixVectorProduct(const matrix_type & A,
			 const vector_type & b,
			 vector_type & c,
			 bool transposeA = false)
{
  assert( A.isFillingCompleted() );
  assert( A.globalCols() == b.globalSize() );
  A.data()->Multiply(transposeA, *b.data(), *c.data());
}

  

  
} // end namespace core
#endif
