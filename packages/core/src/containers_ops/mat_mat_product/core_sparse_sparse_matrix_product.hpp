
#ifndef CORE_MATRIX_OPERATIONS_SPARSE_SPARSE_MATRIX_PRODUCT_HPP_
#define CORE_MATRIX_OPERATIONS_SPARSE_SPARSE_MATRIX_PRODUCT_HPP_

#include "../core_ops_meta.hpp"
#include "../../matrix/core_matrix_meta.hpp"

namespace rompp{ namespace core{ namespace ops{
  
/*---------------------------------------------------------
C = A B
- A is sparse matrix from eigen 
- B is sparse matrix from eigen 
---------------------------------------------------------*/
  
template <typename A_t, typename B_t, typename C_t,
  core::meta::enable_if_t<
    core::meta::is_eigen_sparse_matrix_wrapper<A_t>::value &&
    core::meta::is_eigen_sparse_matrix_wrapper<B_t>::value &&
    core::meta::is_eigen_sparse_matrix_wrapper<C_t>::value && 
    core::meta::wrapper_triplet_have_same_scalar<A_t,B_t,C_t>::value
    > * = nullptr
  >
void product(const A_t & A, const B_t & B, C_t & C){

  assert(C.rows() == A.rows());
  assert(C.cols() == B.cols());
  (*C.data()) = (*A.data()) * (*B.data());
}

template <typename A_t, typename B_t, 
    core::meta::enable_if_t<
      core::meta::is_eigen_sparse_matrix_wrapper<A_t>::value &&
      core::meta::is_eigen_sparse_matrix_wrapper<B_t>::value &&
      core::meta::wrapper_pair_have_same_scalar<A_t,B_t>::value
      > * = nullptr
    >
A_t product(const A_t & A, const B_t & B){
  
  A_t C(A.rows(),B.cols());
  product(A, B, C);
  return C;
}



}}} // end namespace rompp::core::ops
#endif























// /*---------------------------------------------------
//   C = A * B
//   A: epetra sparse matrix
//   B: epetra sparse matrix
//   C: epetra sparse matrix
// ---------------------------------------------------*/

// template <typename mat_type,
// 	  typename std::enable_if<
// 	    details::traits<mat_type>::isEpetra &&
// 	    details::traits<mat_type>::is_sparse
// 	    >::type * = nullptr
// 	  >
// auto product(const mat_type & A,
// 	     const mat_type & B,
// 	     bool transposeA,
// 	     bool transposeB,
// 	     bool call_filingIsCompleted_on_result = true)
// {  
//   assert( A.isFillingCompleted() );
//   assert( B.isFillingCompleted() );

//   // guess number of non zeros
//   auto maxNonzB = B.data()->GlobalMaxNumEntries();

//   // rowmap of C is the range map of A if we use A
//   // rowmap of C is the domain map of A if we use A^T
//   auto & Crowmap = transposeA ? A.getDomainDataMap()
//     : A.getRangeDataMap();
//   mat_type C(Crowmap, maxNonzB);

//   core::mat_ops::product(A, B, C, transposeA, transposeB,
//    		      call_filingIsCompleted_on_result);
//   return C;

// }//end fnc




// /*---------------------------------------------------
//   C = A * B
//   A: epetra CRS matrix
//   B: epetra CRS matrix
//   C: epetra CRS matrix
// ---------------------------------------------------*/

// template <typename mat_type,
// 	  typename std::enable_if<
// 	    details::traits<mat_type>::isEpetra &&
// 	    details::traits<mat_type>::is_sparse
// 	    >::type * = nullptr
// 	  >
// void product(const mat_type & A,
// 	     const mat_type & B,
// 	     mat_type & C,
// 	     bool transposeA,
// 	     bool transposeB,
// 	     bool call_filingIsCompleted_on_result = true){

//   // /* From trilinos docs: 
//   //    Given Epetra_CrsMatrix objects A, B and C, form the product C = A*B.
//   //    In a parallel setting, A and B need not have matching distributions, 
//   //    but C needs to have the same row-map as A.

//   //    Parameters:
//   //    *A Input, must already have had 'FillComplete()' called.
//   //    *transposeA Input, whether to use transpose of matrix A.
//   //    *B Input, must already have had 'FillComplete()' called.
//   //    *transposeB Input, whether to use transpose of matrix B.
//   //    *C Result. On entry to this method, it doesn't matter whether 
//   //    FillComplete() has already been called on C or not. If it has, 
//   //    then C's graph must already contain all nonzero locations that 
//   //    will be produced when forming the product A*B. On exit, C.FillComplete() 
//   //    will have been called, unless the last argument to this function 
//   //    is specified to be false.
//   //    * call_FillComplete_on_result Optional argument, defaults to true. 
//   //    Power users may specify this argument to be false if they *DON'T* 
//   //    want this function to call C.FillComplete. (It is often useful 
//   //    to allow this function to call C.FillComplete, in cases where one 
//   //    or both of the input matrices are rectangular and it is not 
//   //    trivial to know which maps to use for the domain- and range-maps.)
//   //  */
  
//   assert( A.isFillingCompleted() );
//   assert( B.isFillingCompleted() );
//   auto & rangeMapAB = transposeA ? A.getDomainDataMap() : A.getRowDataMap();
//   auto & rowMapC = C.getRowDataMap();
//   assert( rowMapC.SameAs(rangeMapAB) );

//   if ( C.isFillingCompleted() ){
//     assert( C.hasSameRangeDataMapAs(A) );
//     assert( C.hasSameDomainDataMapAs(B) );
//     EpetraExt::MatrixMatrix::Multiply(*A.data(),
// 		   transposeA,
// 		   *B.data(), transposeB,
// 		   *C.data(),
// 		   call_filingIsCompleted_on_result);
//   }
//   else{
//     EpetraExt::MatrixMatrix::Multiply(*A.data(), transposeA,
// 		  *B.data(), transposeB,
// 		  *C.data(),
// 		  call_filingIsCompleted_on_result);
//     assert( C.data()->Filled() );
//   }

// }//end fnc
