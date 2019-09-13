/*
//@HEADER
// ************************************************************************
//
// containers_sparse_matrix_transpose.hpp
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

#ifndef CONTAINERS_MATRIX_OPERATIONS_SPARSE_MATRIX_TRANSPOSE_HPP_
#define CONTAINERS_MATRIX_OPERATIONS_SPARSE_MATRIX_TRANSPOSE_HPP_

#include "../../meta/containers_matrix_meta.hpp"
#include "../concrete/containers_matrix_sparse_sharedmem_eigen.hpp"
#include "EpetraExt_Transpose_RowMatrix.h"
//#include <Epetra_RowMatrixTransposer.h>

namespace pressio{
namespace containers{
namespace mat_ops{
  
/*-----------------------------------------------------
  EPETRA CRSMATRIX
----------------------------------------------------- */
template <typename mat_type,
	  ::pressio::mpl::enable_if_t<
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
  // containers::Matrix<nat_t> res( *transA );
  // assert( res.isFillingCompleted() );
  // delete transA;

  //-----------
  // method 2
  //-----------
  EpetraExt::RowMatrix_Transpose transposer;
  Epetra_CrsMatrix & transA =
    dynamic_cast<Epetra_CrsMatrix &>(transposer(*A.data()));
  containers::Matrix<nat_t> res( transA );
  assert( res.isFillingCompleted() );
  return res;
}

  
/*-----------------------------------------------------
  EIGEN SPARSE
----------------------------------------------------- */
template <typename mat_type,
	  ::pressio::mpl::enable_if_t<
	    details::traits<mat_type>::isEigen &&
	    details::traits<mat_type>::is_sparse 
	    > * = nullptr>
auto transpose(const mat_type & A)
{
  mat_type res(A.data()->transpose());
  return res;
}


  
} // end namespace mat_ops  
} // end namespace containers
}//end namespace pressio
#endif
