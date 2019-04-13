
#ifndef CORE_MATRIX_OPERATIONS_SPARSE_MATRIX_TRANSPOSE_HPP_
#define CORE_MATRIX_OPERATIONS_SPARSE_MATRIX_TRANSPOSE_HPP_

#include "../../meta/core_matrix_meta.hpp"
#include "../concrete/core_matrix_sparse_sharedmem_eigen.hpp"
#include "EpetraExt_Transpose_RowMatrix.h"
//#include <Epetra_RowMatrixTransposer.h>

namespace rompp{
namespace core{
namespace mat_ops{
  
/*-----------------------------------------------------
  EPETRA CRSMATRIX
----------------------------------------------------- */
template <typename mat_type,
	  ::rompp::mpl::enable_if_t<
	    details::traits<mat_type>::isEpetra &&
	    details::traits<mat_type>::is_sparse 
	    > * = nullptr>
auto transpose(mat_type & A
	       /*nonconst & because it is needed 
		 by rowmatrix transposer*/) 
{
  using nat_t = typename details::traits<mat_type>::wrapped_t;

  //-----------
  // method 1
  //-----------
  // Epetra_RowMatrixTransposer transposer(A.data());
  // nat_t * transA;
  // transposer.CreateTranspose(false, transA);
  // core::Matrix<nat_t> res( *transA );
  // assert( res.isFillingCompleted() );
  // delete transA;

  //-----------
  // method 2
  //-----------
  EpetraExt::RowMatrix_Transpose transposer;
  Epetra_CrsMatrix & transA =
    dynamic_cast<Epetra_CrsMatrix &>(transposer(*A.data()));
  core::Matrix<nat_t> res( transA );
  assert( res.isFillingCompleted() );
  return res;
}

  
/*-----------------------------------------------------
  EIGEN SPARSE
----------------------------------------------------- */
template <typename mat_type,
	  ::rompp::mpl::enable_if_t<
	    details::traits<mat_type>::isEigen &&
	    details::traits<mat_type>::is_sparse 
	    > * = nullptr>
auto transpose(const mat_type & A)
{
  mat_type res(A.data()->transpose());
  return res;
}


  
} // end namespace mat_ops  
} // end namespace core
}//end namespace rompp
#endif
