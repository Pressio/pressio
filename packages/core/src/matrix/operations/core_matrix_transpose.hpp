
#ifndef CORE_MATRIX_OPERATIONS_MATRIX_TRANSPOSE_HPP_
#define CORE_MATRIX_OPERATIONS_MATRIX_TRANSPOSE_HPP_

#include "../../meta/core_matrix_meta.hpp"
#include "../concrete/core_matrix_dense_sharedmem_eigen.hpp"
#include "../concrete/core_matrix_sparse_sharedmem_eigen.hpp"
#include "EpetraExt_Transpose_RowMatrix.h"
//#include <Epetra_RowMatrixTransposer.h>

namespace core{


/*=====================================
This is a list supported:

Legend: 
epCRS = epetra CRS matrix
epDM = epetra dense matrix

* transpose(epCRS) (ok)
* transpose(eigenDense) (ok)
* transpose(eigenSparse) (ok)

* transpose(epDM) (missing)

===================================*/


/*-----------------------------------------------------
  EPETRA DENSE MATRIX
----------------------------------------------------- */
template <typename mat_type,
	  typename std::enable_if<
	    details::traits<mat_type>::isEpetra &&
	    details::traits<mat_type>::isDense
	    >::type * = nullptr>
auto transpose(const mat_type & A) 
{
  // auto & mapA = A.getDataMap();
  auto const & mpicomm = A.commCRef();
    
  // let's ensure that the transposed matrix cna be distributed
  // over the current processes, so the nuber of of cols of A
  // needs to be greater than the number of MPI processes running 
  const int numProc = mpicomm.NumProc();
  if( numProc > A.globalRows() ){
    std::cout << " ERROR: trying to transpose a dense dist matrix \
with number of columns too small to be distributed over current communicator \n";}
  assert( numProc < A.globalRows() );

  Epetra_Map mapT(A.globalCols(), 0, mpicomm);  
  mat_type C( mapT, A.globalRows() );
  C.setZero();

  // missing impl
  
  return C;
}


/*-----------------------------------------------------
  EPETRA CRSMATRIX
----------------------------------------------------- */
template <typename mat_type,
	  typename std::enable_if<
	    details::traits<mat_type>::isEpetra &&
	    details::traits<mat_type>::isSparse 
	    >::type * = nullptr>
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


