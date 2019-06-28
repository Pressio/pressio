
#ifndef ALGEBRA_MATRIX_DENSE_TO_SPARSE_HPP_
#define ALGEBRA_MATRIX_DENSE_TO_SPARSE_HPP_

#include "../../meta/algebra_matrix_meta.hpp"
#include "../concrete/algebra_matrix_dense_sharedmem_eigen.hpp"
#include "../concrete/algebra_matrix_sparse_sharedmem_eigen.hpp"
#include "../concrete/algebra_matrix_dense_distributed_epetra.hpp"
#include "../concrete/algebra_matrix_sparse_distributed_epetra.hpp"

/*=====================================
Transform a dense matrix to sparse 
===================================*/


namespace rompp{
namespace algebra{
  
/*---------------------------------------------------
A => B where Epetra dense A to CRS B.
Here we pass the matrix and domain and range maps.
---------------------------------------------------*/
template <typename mat_type,
	  typename std::enable_if<
	    details::traits<mat_type>::isEpetra &&
	    details::traits<mat_type>::is_dense
	    >::type * = nullptr
	  >
auto denseToSparse(const mat_type & A,
		   const Epetra_Map & domain_map,
		   const Epetra_Map & range_map)
{
  
  // non zeros per row
  const auto nzPerRow = A.globalCols();
  // row map: rememeber that A is a dense matrix
  auto & rowmap = static_cast<const Epetra_Map &>(A.getDataMap());

  // now create the target CRS matrix
  using res_t = algebra::Matrix<Epetra_CrsMatrix>;
  using traits_t = details::traits<res_t>;
  using sc_t = typename traits_t::scalar_t;
  using LO_t = typename traits_t::local_ordinal_t;
  using GO_t = typename traits_t::global_ordinal_t;

  // the CRS matrix
  res_t B(rowmap, nzPerRow);

  // fill it
  const LO_t numMyEl = rowmap.NumMyElements();
  std::vector<sc_t> val(nzPerRow);
  std::vector<GO_t> ind(nzPerRow);
  for (GO_t j=0; j<nzPerRow; j++)
    ind[j]=j;

  std::vector<GO_t> mygel(numMyEl);
  rowmap.MyGlobalElements(mygel.data());  
  for (LO_t i=0; i<numMyEl; i++){
    auto grow = mygel[i];
    for (GO_t j=0; j<nzPerRow; j++){
      val[j]=A(i,j);
    }
    B.insertGlobalValues(grow, nzPerRow, val.data(), ind.data());
  }

  //fill complete
  B.fillingIsCompleted(domain_map, range_map);
  
  return B;
}

  
/*
A => B where Epetra dense A to CRS B.

if just the matrix is passed, then fill complete is 
called on the matrix using the rowmap for both 
domain and range maps. This is what the regular 
fillcomplete() does within trilinos.
*/
template <typename mat_type,
	  typename std::enable_if<
	    details::traits<mat_type>::isEpetra &&
	    details::traits<mat_type>::is_dense
	    >::type * = nullptr
	  >
auto denseToSparse(const mat_type & A)
{

  auto & rowmap = static_cast<const Epetra_Map &>(A.getDataMap());
  Epetra_Map colmap( A.globalCols(), 0, A.commCRef() );
  return denseToSparse(A, colmap, rowmap);
  
}

  
} // end namespace algebra
}//end namespace rompp
#endif
