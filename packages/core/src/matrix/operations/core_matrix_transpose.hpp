
#ifndef CORE_MATRIX_OPERATIONS_MATRIX_TRANSPOSE_HPP_
#define CORE_MATRIX_OPERATIONS_MATRIX_TRANSPOSE_HPP_

#include "../../meta/core_matrix_meta.hpp"
#include "../concrete/core_matrix_dense_serial_eigen.hpp"
#include "../concrete/core_matrix_sparse_serial_eigen.hpp"
#include "EpetraExt_Transpose_RowMatrix.h"
//#include <Epetra_RowMatrixTransposer.h>

namespace core{

/*-----------------------------------------------------
  EPETRA CRSMATRIX
----------------------------------------------------- */
template <typename mat_type,
	  typename std::enable_if<
	    details::traits<mat_type>::isEpetra &&
	    details::traits<mat_type>::isSparse 
	    >::type * = nullptr>
auto transpose(mat_type & A /*nonconst & because it is needed by rowmatrix transposer*/) 
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
  Epetra_CrsMatrix & transA = dynamic_cast<Epetra_CrsMatrix&>(transposer(*A.data()));
  core::Matrix<nat_t> res( transA );
  assert( res.isFillingCompleted() );
  return res;
}

  
/*-----------------------------------------------------
  EIGEN SPARSE
----------------------------------------------------- */
template <typename mat_type,
	  typename std::enable_if<
	    details::traits<mat_type>::isEigen &&
	    details::traits<mat_type>::isSparse 
	    >::type * = nullptr>
auto transpose(const mat_type & A)
{
  mat_type res(A.data()->transpose());
  return res;
}
  
/*-----------------------------------------------------
  EIGEN DENSE DYNAMIC
----------------------------------------------------- */
template <typename mat_type,
	  typename std::enable_if<
	    details::traits<mat_type>::isEigen &&
	    details::traits<mat_type>::isDense &&
	    details::traits<mat_type>::isStatic==0
	    >::type * = nullptr>
mat_type transpose(const mat_type & A)
{
  mat_type res(A);
  res.data()->transposeInPlace();
  return res;
}

/*-----------------------------------------------------
  EIGEN DENSE STATIC
----------------------------------------------------- */
template <typename mat_type,
	  typename std::enable_if<
	    details::traits<mat_type>::isEigen &&
	    details::traits<mat_type>::isDense &&
	    details::traits<mat_type>::isStatic
	    >::type * = nullptr>
auto transpose(const mat_type & A)
{
  using scalar_t = typename details::traits<mat_type>::scalar_t;
  static constexpr int rowsIn
    = details::traits<mat_type>::wrapped_t::RowsAtCompileTime;
  static constexpr int colsIn
    = details::traits<mat_type>::wrapped_t::ColsAtCompileTime;

  using res_t = core::Matrix<Eigen::Matrix<scalar_t, colsIn, rowsIn>>;
  res_t res;
  *res.data() = A.data()->transpose().eval();
  return res;
}
  
  
} // end namespace core
#endif


