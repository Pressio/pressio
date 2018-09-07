
#ifndef CORE_MATRIX_OPERATIONS_SPARSE_SPARSE_MATRIX_PRODUCT_HPP_
#define CORE_MATRIX_OPERATIONS_SPARSE_SPARSE_MATRIX_PRODUCT_HPP_

#include "../../meta/core_matrix_meta.hpp"
#include "../concrete/core_matrix_dense_sharedmem_eigen.hpp"
#include "../concrete/core_matrix_sparse_sharedmem_eigen.hpp"
#include <EpetraExt_MatrixMatrix.h>
#include "TpetraExt_MatrixMatrix.hpp"
#include "Epetra_LocalMap.h"
#include <Epetra_Export.h>
#include <Epetra_Import.h>


namespace core{
namespace mat_ops{

  
/*---------------------------------------------------
-----------------------------------------------------
  C = A * B
  A: epetra CRS matrix
  B: epetra CRS matrix
  C: epetra CRS matrix
-----------------------------------------------------
---------------------------------------------------*/

template <typename mat_type,
	  typename std::enable_if<
	    details::traits<mat_type>::isEpetra &&
	    details::traits<mat_type>::isSparse
	    >::type * = nullptr
	  >
void product(const mat_type & A,
	     const mat_type & B,
	     mat_type & C,
	     bool transposeA,
	     bool transposeB,
	     bool call_filingIsCompleted_on_result = true)
{

  /* From trilinos docs: 

     Given Epetra_CrsMatrix objects A, B and C, form the product C = A*B.
     In a parallel setting, A and B need not have matching distributions, 
     but C needs to have the same row-map as A.

     Parameters:
     *A Input, must already have had 'FillComplete()' called.
     *transposeA Input, whether to use transpose of matrix A.
     *B Input, must already have had 'FillComplete()' called.
     *transposeB Input, whether to use transpose of matrix B.
     *C Result. On entry to this method, it doesn't matter whether 
     FillComplete() has already been called on C or not. If it has, 
     then C's graph must already contain all nonzero locations that 
     will be produced when forming the product A*B. On exit, C.FillComplete() 
     will have been called, unless the last argument to this function 
     is specified to be false.
     * call_FillComplete_on_result Optional argument, defaults to true. 
     Power users may specify this argument to be false if they *DON'T* 
     want this function to call C.FillComplete. (It is often useful 
     to allow this function to call C.FillComplete, in cases where one 
     or both of the input matrices are rectangular and it is not 
     trivial to know which maps to use for the domain- and range-maps.)
   */
  
  assert( A.isFillingCompleted() );
  assert( B.isFillingCompleted() );

  auto & rangeMapAB = transposeA ? A.getDomainDataMap() : A.getRowDataMap();
  auto & rowMapC = C.getRowDataMap();
  assert( rowMapC.SameAs(rangeMapAB) );

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
    assert( C.data()->Filled() );
  }

}//end fnc


/*---------------------------------------------------
-----------------------------------------------------
  C = A * B
  A: epetra sparse matrix
  B: epetra sparse matrix
  C: epetra sparse matrix
-----------------------------------------------------
---------------------------------------------------*/

template <typename mat_type,
	  typename std::enable_if<
	    details::traits<mat_type>::isEpetra &&
	    details::traits<mat_type>::isSparse
	    >::type * = nullptr
	  >
auto product(const mat_type & A,
	     const mat_type & B,
	     bool transposeA,
	     bool transposeB,
	     bool call_filingIsCompleted_on_result = true)
{  
  assert( A.isFillingCompleted() );
  assert( B.isFillingCompleted() );

  // guess number of non zeros
  auto maxNonzB = B.data()->GlobalMaxNumEntries();

  // rowmap of C is the range map of A if we use A
  // rowmap of C is the domain map of A if we use A^T
  auto & Crowmap = transposeA ? A.getDomainDataMap() : A.getRangeDataMap();
  mat_type C(Crowmap, maxNonzB);

  core::mat_ops::product(A, B, C, transposeA, transposeB,
   		      call_filingIsCompleted_on_result);
  return C;

}//end fnc


} // end namespace mat_ops
} // end namespace core
#endif
