
#ifndef CORE_MATRIX_OPERATIONS_SPARSE_MATRIX_VECTOR_PRODUCT_HPP_
#define CORE_MATRIX_OPERATIONS_SPARSE_MATRIX_VECTOR_PRODUCT_HPP_

#include "core_ops_meta.hpp"
#include "../vector/core_vector_meta.hpp"
#include "../matrix/core_matrix_meta.hpp"

namespace rompp{
namespace core{
namespace ops{
    
/*---------------------------------------------------------
c = A b
- A is matrix from eigen
- b is vector from eigen
- c is an eigen vector storing the result 
---------------------------------------------------------*/
  
template <typename mat_t, typename vec1_t, typename vec2_t,
	  core::meta::enable_if_t<
	    core::meta::is_eigen_sparse_matrix_wrapper<mat_t>::value &&
	    core::meta::is_eigen_vector_wrapper<vec1_t>::value &&
	    core::meta::is_eigen_vector_wrapper<vec2_t>::value
	    > * = nullptr
	  >
void product(const mat_t & A, const vec1_t & b, vec2_t & c){

  assert(A.cols() == b.size());
  assert(c.size() == A.rows());
  (*c.data()) = (*A.data()) * (*b.data());
}


  
} // end namespace ops
} // end namespace core
}//end namespace rompp
#endif




// /*---------------------------------------------------------
//   EPETRA 
//   c = A b , 
//   - A = crs matrix 
//   - b = SINGLE vector
//   - c is epetra vector
// -----------------------------------------------------------*/

// template <typename matrix_type,
// 	  typename vector_type,
// 	  core::meta::enable_if_t<
// 	    details::traits<matrix_type>::is_matrix==1 &&
// 	    details::traits<matrix_type>::isEpetra==1 &&
// 	    details::traits<matrix_type>::is_sparse==1 &&
// 	    details::traits<vector_type>::is_vector==1 &&
// 	    details::traits<vector_type>::isEpetra==1
// 	    > * = nullptr
// 	  >
// auto product(const matrix_type & A,
// 	     const vector_type & b,
// 	     bool transposeA = false){

//   assert( A.isFillingCompleted() );
//   assert( A.globalCols() == b.globalSize() );
//   vector_type c( A.getRangeDataMap() );
//   A.data()->Multiply(transposeA, *b.data(), *c.data());
//   return c;
// }


// template <typename matrix_type,
// 	  typename vector_type,
// 	  core::meta::enable_if_t<
// 	    details::traits<matrix_type>::is_matrix==1 &&
// 	    details::traits<matrix_type>::isEpetra==1 &&
// 	    details::traits<matrix_type>::is_sparse==1 &&
// 	    details::traits<vector_type>::is_vector==1 &&
// 	    details::traits<vector_type>::isEpetra==1
// 	    > * = nullptr
// 	  >
// void product(const matrix_type & A,
// 	     const vector_type & b,
// 	     vector_type & c,
// 	     bool transposeA = false){

//   assert( A.isFillingCompleted() );
//   assert( A.globalCols() == b.globalSize() );
//   A.data()->Multiply(transposeA, *b.data(), *c.data());
// }
