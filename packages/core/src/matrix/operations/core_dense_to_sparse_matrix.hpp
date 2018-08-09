
#ifndef CORE_MATRIX_DENSE_TO_SPARSE_HPP_
#define CORE_MATRIX_DENSE_TO_SPARSE_HPP_

#include "../../meta/core_matrix_meta.hpp"
#include "../concrete/core_matrix_dense_serial_eigen.hpp"
#include "../concrete/core_matrix_sparse_serial_eigen.hpp"
#include "../concrete/core_matrix_dense_distributed_epetra.hpp"
#include "../concrete/core_matrix_sparse_distributed_epetra.hpp"


/*=====================================
Transform a dense matrix to sparse 
===================================*/

namespace core{
  
/*---------------------------------------------------
-----------------------------------------------------
A => B
Epetra dense A to CRS B
-----------------------------------------------------
---------------------------------------------------*/
template <typename mat_type,
	  typename std::enable_if<
	    details::traits<mat_type>::isEpetra &&
	    details::traits<mat_type>::isDense
	    >::type * = nullptr
	  >
auto denseToSparse(const mat_type & A)
{
  // non zeros per row, single value
  const auto nzPerRow = A.globalCols();
  // row map
  auto & dmap = A.getDataMap();
 
  // target CRS matrix
  using res_t = core::Matrix<Epetra_CrsMatrix>;
  using traits_t = details::traits<res_t>;
  using sc_t = typename traits_t::scalar_t;
  using LO_t = typename traits_t::local_ordinal_t;
  using GO_t = typename traits_t::global_ordinal_t;

  res_t B(dmap, nzPerRow);

  const LO_t numMyEl = dmap.NumMyElements();
  std::vector<sc_t> val(nzPerRow);
  std::vector<GO_t> ind(nzPerRow);
  for (GO_t j=0; j<nzPerRow; j++)
    ind[j]=j;

  std::vector<GO_t> mygel(numMyEl);
  dmap.MyGlobalElements(mygel.data());
  

  for (LO_t i=0; i<numMyEl; i++)
  {
    auto grow = mygel[i];
    for (GO_t j=0; j<nzPerRow; j++){
      val[j]=A(i,j);
    }
    B.insertGlobalValues(grow, nzPerRow, val.data(), ind.data());
  }  
  return B;
}


   
} // end namespace core
#endif
