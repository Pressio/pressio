
#ifndef CORE_MATRIX_OPERATIONS_SPARSE_DENSE_MATRIX_PRODUCT_HPP_
#define CORE_MATRIX_OPERATIONS_SPARSE_DENSE_MATRIX_PRODUCT_HPP_

#include "../../meta/core_matrix_meta.hpp"
#include "../concrete/core_matrix_dense_sharedmem_eigen.hpp"
#include "../concrete/core_matrix_sparse_sharedmem_eigen.hpp"
#include <EpetraExt_MatrixMatrix.h>
#include "TpetraExt_MatrixMatrix.hpp"
#include "Epetra_LocalMap.h"
#include <Epetra_Export.h>
#include <Epetra_Import.h>

namespace rompp{
namespace core{
namespace mat_ops{

  
// /*-----------------------------------------------------
// -------------------------------------------------------
//   C = A * B

//   A: epetra CRS matrix wrapper
//   B: epetra dense matrix wrapper
//   C: epetra dense matrix wrapper
// -------------------------------------------------------
// -----------------------------------------------------*/
// template <typename mat_sp_type,
// 	  typename mat_ds_type,
// 	  typename std::enable_if<
// 	    details::traits<mat_ds_type>::isEpetra &&
// 	    details::traits<mat_ds_type>::is_dense &&
// 	    details::traits<mat_sp_type>::isEpetra &&
// 	    details::traits<mat_sp_type>::is_sparse
// 	    >::type * = nullptr
// 	  >
// auto product(const mat_sp_type & A,
// 			 const mat_ds_type & B)
// {
//   // get row map of A
//   auto & mapA = A.getRangeDataMap();
//   mat_ds_type C( mapA, B.globalCols() );
//   C.setZero();
//   A.data()->Multiply(false, *B.data(), *C.data());
//   return C;
// }


  
} // end namespace mat_ops
} // end namespace core
}//end namespace rompp
#endif
