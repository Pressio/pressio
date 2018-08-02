
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

template <typename mat_type,
	  typename std::enable_if<
	    details::traits<mat_type>::isEpetra &&
	    details::traits<mat_type>::isSparse
	    >::type * = nullptr
	  >
void matrixMatrixProduct(const mat_type & A,
			 const mat_type & B,
			 mat_type & C,
			 bool transposeA,
			 bool transposeB,
			 bool call_filingIsCompleted_on_result = true)
{

  assert( A.isFillingCompleted() );
  assert( B.isFillingCompleted() );
  assert( C.hasSameRowDataMapAs(A) );

  if ( C.isFillingCompleted() ){
    assert( C.hasSameRangeDataMapAs(A) );
    assert( C.hasSameDomainDataMapAs(B) );
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

  

/*-----------------------------------------------------
-------------------------------------------------------
  C = A * B

  A: epetra dense matrix
  B: epetra dense matrix
  C: epetra dense matrix
-------------------------------------------------------
-----------------------------------------------------*/
template <typename mat_type,
	  typename std::enable_if<
	    details::traits<mat_type>::isEpetra &&
	    details::traits<mat_type>::isDense
	    >::type * = nullptr
	  >
auto matrixMatrixProduct(const mat_type & A,
			 const mat_type & B)
{
  assert( A.globalCols() == B.globalRows() );

  /* I tried here to use the Multiply method of MultiVectors but it 
     does not seem to work as expected. When A,B,C are all distributed, I don't get 
     the right result. So we need to figure out why. 
     Nonetheless, for now put a placeholder doing nothing. Just create 
     the result dense matrix and then fill it with zeros. 
     Maybe we should just put here for now a simple implementation.
  */
  
  // get row map of A
  auto & mapA = A.getDataMap();
  mat_type C( mapA, B.globalCols() );
  C.setZero();
  return C;
}



/*-----------------------------------------------------
-------------------------------------------------------
  C = A * B

  A: epetra CRS matrix wrapper
  B: epetra dense matrix wrapper
  C: epetra dense matrix wrapper
-------------------------------------------------------
-----------------------------------------------------*/
template <typename mat_sp_type,
	  typename mat_ds_type,
	  typename std::enable_if<
	    details::traits<mat_ds_type>::isEpetra &&
	    details::traits<mat_ds_type>::isDense &&
	    details::traits<mat_sp_type>::isEpetra &&
	    details::traits<mat_sp_type>::isSparse
	    >::type * = nullptr
	  >
auto matrixMatrixProduct(const mat_sp_type & A,
			 const mat_ds_type & B)
{
  assert( A.globalCols() == B.globalRows() );

  std::cout << " OKKKK \n"; 

  // get row map of A
  auto & mapA = A.getRangeDataMap();
  mat_ds_type C( mapA, B.globalCols() );
  C.setZero();

  A.data()->Multiply(false, *B.data(), *C.data());
  return C;
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
