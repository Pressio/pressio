
#ifndef CORE_MATRIX_OPERATIONS_MATRIX_MATRIX_PRODUCT_HPP_
#define CORE_MATRIX_OPERATIONS_MATRIX_MATRIX_PRODUCT_HPP_

#include "../../meta/core_matrix_meta.hpp"
#include "../concrete/core_matrix_dense_serial_eigen.hpp"
#include "../concrete/core_matrix_sparse_serial_eigen.hpp"
#include <EpetraExt_MatrixMatrix.h>
#include "TpetraExt_MatrixMatrix.hpp"

namespace core{

/*-----------------------------------------------------
  C = A * B
   A: epetra sparse matrix
   B: epetra sparse matrix
   C: epetra sparse matrix
----------------------------------------------------- */
template <typename mat_type>
void
matrixMatrixProduct(const mat_type & A,
		    const mat_type & B,
		    mat_type & C,
		    bool transposeA = false,
		    bool transposeB = false,
		    bool call_filingIsCompleted_on_result = true,
		    typename std::enable_if<
		     details::traits<mat_type>::isEpetra &&
		     details::traits<mat_type>::isSparse
		    >::type * = nullptr)
{

  assert( A.isFillingCompleted() );
  assert( B.isFillingCompleted() );
  assert( C.rowMap() == A.rowMap() );

  if ( C.isFillingCompleted() ){
    assert( C.rangeMap() == A.rangeMap() );
    assert( C.domainMap() == B.domainMap() );
    EpetraExt::MatrixMatrix::Multiply(*A.data(), transposeA,
				      *B.data(), transposeB,
				      *C.data(),
				      call_filingIsCompleted_on_result);
  }
  else{
    EpetraExt::MatrixMatrix::Multiply(*A.data(), transposeA,
				      *B.data(), transposeB,
				      *C.data(),
				      call_filingIsCompleted_on_result);
    C.fillingIsCompleted();
  }
}


  



// /*-----------------------------------------------------
//    TPETRA
//   C = A * B
//    A: tpetra sparse matrix
//    B: tpetra sparse matrix
//    C: tpetra sparse matrix
// ----------------------------------------------------- */
// template <typename mat_type>
// void
// matrixMatrixProduct(const mat_type & A,
// 		    const mat_type & B,
// 		    mat_type & C,
// 		    bool transposeA = false,
// 		    bool transposeB = false,
// 		    bool call_filingIsCompleted_on_result = true,
// 		    typename std::enable_if<
// 		    meta::is_matrix_sparse_distributed_tpetra<mat_type>::value
// 		    >::type * = nullptr)
// {
//   assert( A.isFillingCompleted() );
//   assert( B.isFillingCompleted() );
//   assert( C.rowMap() == A.rowMap() );
//   assert( A.globalCols() == B.globalRows() );

//   if ( C.isFillingCompleted() ){
//     assert( C.rangeMap() == A.rangeMap() );
//     assert( C.domainMap() == B.domainMap() );
//     EpetraExt::MatrixMatrix::Multiply(*A.data(), transposeA,
// 				      *B.data(), transposeB,
// 				      *C.data(),
// 				      call_filingIsCompleted_on_result);
//   }
//   else{
//     EpetraExt::MatrixMatrix::Multiply(*A.data(), transposeA,
// 				      *B.data(), transposeB,
// 				      *C.data(),
// 				      call_filingIsCompleted_on_result);
//     C.fillingIsCompleted();
//   }
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
