
#ifndef CORE_MATRIX_OPERATIONS_MATRIX_MATRIX_PRODUCT_HPP_
#define CORE_MATRIX_OPERATIONS_MATRIX_MATRIX_PRODUCT_HPP_

#include "../../meta/core_matrix_meta.hpp"
#include "../concrete/core_matrix_dense_sharedmem_eigen.hpp"
#include "../concrete/core_matrix_sparse_sharedmem_eigen.hpp"
#include <EpetraExt_MatrixMatrix.h>
#include "TpetraExt_MatrixMatrix.hpp"
#include "Epetra_LocalMap.h"
#include <Epetra_Export.h>
#include <Epetra_Import.h>


/*=====================================
This is a list of products supported:

Legend: 
epCRS = epetra CRS matrix
epDM = epetra dense matrix

* epCRS = epCRS * epCRS (ok)
* epDM = epCRS * epDM (ok)
* epDM = epDM * epDM (done but not optimal)

* epDM = epDM * epCRS (missing)

===================================*/


namespace core{

  
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
void matrixMatrixProduct(const mat_type & A,
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
auto matrixMatrixProduct(const mat_type & A,
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

  matrixMatrixProduct(A, B, C, transposeA, transposeB,
   		      call_filingIsCompleted_on_result);
  return C;

}//end fnc



  
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

  const auto BGSize = B.globalRows();
  assert( A.globalCols() == BGSize );
  
  /* I tried here to use the Multiply method of MultiVectors 
     but it does not seem to work as expected. 
     When A,B are all distributed, I don't get 
     the right result. So we need to figure out why. 

     Only solution that worked is to do this trick: 
        B is distributed -> import into B replicated -> do multiply

     We should find out why not working for fully distributed case
  */  

  // define local map
  Epetra_LocalMap locMap( BGSize, 0, B.commCRef() );
  // define replicated B
  Epetra_MultiVector BRep(locMap, B.globalCols());
    
  // get distributed map
  auto & srcMap = B.getDataMap();
  // define importer: Epetra_Import(targetMap, sourceMap)
  Epetra_Import globToLocalImporter(locMap, srcMap);

  // import global -> local
  BRep.Import(*B.data(), globToLocalImporter, Insert);
  
  mat_type C( A.getDataMap(), B.globalCols() );
  C.data()->Multiply( 'N','N', 1.0,  *A.data(), BRep, 0.0 );  
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
  // get row map of A
  auto & mapA = A.getRangeDataMap();
  mat_ds_type C( mapA, B.globalCols() );
  C.setZero();
  A.data()->Multiply(false, *B.data(), *C.data());
  return C;
}
  
  
} // end namespace core
#endif
