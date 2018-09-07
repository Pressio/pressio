
#ifndef CORE_MATRIX_OPERATIONS_DENSE_DENSE_MATRIX_PRODUCT_HPP_
#define CORE_MATRIX_OPERATIONS_DENSE_DENSE_MATRIX_PRODUCT_HPP_

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
auto product(const mat_type & A,
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
  
} // end namespace mat_ops
} // end namespace core
#endif
