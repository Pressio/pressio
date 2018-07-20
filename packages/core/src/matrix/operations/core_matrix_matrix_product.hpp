
#ifndef CORE_MATRIX_MATRIX_MATRIX_PRODUCT_HPP_
#define CORE_MATRIX_MATRIX_MATRIX_PRODUCT_HPP_

#include "../../matrix/meta/core_matrix_meta.hpp"
#include "../concrete/core_matrix_dense_serial_eigen.hpp"
#include "../concrete/core_matrix_sparse_serial_eigen.hpp"

namespace core{

// /*-----------------------------------------------------
//   C = A * B
// ----------------------------------------------------- */
// template <typename mat_type,
// 	  typename mat_type_b,
// 	  typename result_type>
// void matrixMatrixProduct(const mat_type & A,
// 			 const mat_type_b & B,
// 			 result_type & C,
// 			 typename std::enable_if<
// 			   // A must be a matrix and from eigen
// 			   details::traits<mat_type>::isMatrix &&
// 			   details::traits<mat_type>::isEigen &&
// 			   //A must be a matrix and from eigen
// 			   details::traits<mat_type_b>::isMatrix &&
// 			   details::traits<mat_type_b>::isEigen &&
// 			   // res must be from eigen
// 			   details::traits<result_type>::isMatrix &&
// 			   details::traits<result_type>::isEigen &&
// 			   // we need to have matching scalar types
// 			   std::is_same<
// 			 typename details::traits<mat_type>::scalar_t,
// 			 typename details::traits<mat_type_b>::scalar_t
// 			   >::value &&
// 			   std::is_same<
// 			 typename details::traits<mat_type>::scalar_t,
// 			 typename details::traits<result_type>::scalar_t
// 			   >::value 			 
// 			 >::type * = nullptr)
// {
//   assert(A.cols() == B.rows());
//   assert(C.rows() == A.rows());
//   assert(C.cols() == B.cols());

//   (*C.data()) = (*A.data()) * (*B.data());
// }


// /*-----------------------------------------------------
//   C = A * B * C
// ----------------------------------------------------- */
// template <typename mat_type_a,
// 	  typename mat_type_b,
// 	  typename result_type>
// void matrixMatrixProduct(const mat_type_a & A,
// 			 const mat_type_b & B,
// 			 result_type & C,
// 			 typename std::enable_if<
// 			   // A must be a matrix and from eigen
// 			   details::traits<mat_type_a>::isMatrix &&
// 			   details::traits<mat_type_a>::isEigen &&
// 			   //A must be a matrix and from eigen
// 			   details::traits<mat_type_b>::isMatrix &&
// 			   details::traits<mat_type_b>::isEigen &&
// 			   // res must be from eigen
// 			   details::traits<result_type>::isMatrix &&
// 			   details::traits<result_type>::isEigen &&
// 			   // we need to have matching scalar types
// 			   std::is_same<
// 			 typename details::traits<mat_type_a>::scalar_t,
// 			 typename details::traits<mat_type_b>::scalar_t
// 			   >::value &&
// 			   std::is_same<
// 			 typename details::traits<mat_type_a>::scalar_t,
// 			 typename details::traits<result_type>::scalar_t
// 			   >::value 			 
// 			 >::type * = nullptr)
// {
//   assert(A.cols() == B.rows());
//   assert(C.rows() == A.rows());
//   assert(C.cols() == B.cols());

//   (*C.data()) = (*A.data()) * (*B.data());
// }


  
} // end namespace core
#endif
