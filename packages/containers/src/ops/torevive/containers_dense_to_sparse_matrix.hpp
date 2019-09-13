/*
//@HEADER
// ************************************************************************
//
// containers_dense_to_sparse_matrix.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef CONTAINERS_MATRIX_DENSE_TO_SPARSE_HPP_
#define CONTAINERS_MATRIX_DENSE_TO_SPARSE_HPP_

#include "../../meta/containers_matrix_meta.hpp"
#include "../concrete/containers_matrix_dense_sharedmem_eigen.hpp"
#include "../concrete/containers_matrix_sparse_sharedmem_eigen.hpp"
#include "../concrete/containers_matrix_dense_distributed_epetra.hpp"
#include "../concrete/containers_matrix_sparse_distributed_epetra.hpp"

/*=====================================
Transform a dense matrix to sparse 
===================================*/


namespace pressio{
namespace containers{
  
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
  using res_t = containers::Matrix<Epetra_CrsMatrix>;
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

  
} // end namespace containers
}//end namespace pressio
#endif
