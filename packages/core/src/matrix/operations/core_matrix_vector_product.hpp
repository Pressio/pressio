
#ifndef CORE_MATRIX_MATRIX_VECTOR_PRODUCT_HPP_
#define CORE_MATRIX_MATRIX_VECTOR_PRODUCT_HPP_

#include "../../vector/meta/core_vector_meta.hpp"
#include "../../vector/concrete/core_vector_serial_eigen.hpp"
#include "../../matrix/meta/core_matrix_meta.hpp"
#include "../concrete/core_matrix_dense_serial_eigen.hpp"
#include "../concrete/core_matrix_sparse_serial_eigen.hpp"

namespace core{

//-----------------------------------------------------
// do: res = A b , where A is matrix and b is vector
// this first case is for when the vector is   
//-----------------------------------------------------  
template <typename matrix_type,
	  typename vector_type,
	  typename result_type>
void matrixVectorProduct(const matrix_type & A,
			 const vector_type & b,
			 result_type & res,
			 typename std::enable_if<
			   // A must be a matrix and from eigen
			   details::traits<matrix_type>::isMatrix &&
			   details::traits<matrix_type>::isEigen &&
			   // b must be a vector and from eigen
			   details::traits<vector_type>::isVector &&
			   details::traits<vector_type>::isEigen &&
			   // res must be a vector and from eigen
			   details::traits<result_type>::isVector &&
			   details::traits<result_type>::isEigen &&
			   // we need to have matching scalar types
			   std::is_same<typename
			                details::traits<matrix_type>::scalar_t,
			                typename
			                details::traits<vector_type>::scalar_t
			               >::value 
			 // // if the vector b is a static vector, need to make sure
			 // // that result has compatible size too
			 // std::conditional<details::traits<vector_type>::isStatic &&
			 //                  >::type
			 >::type * = nullptr)
{
  assert(A.cols() == b.size());
  assert(res.size() == A.rows());
  (*res.data()) = (*A.data()) * (*b.data());
}

} // end namespace core
#endif
