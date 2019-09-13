/*
//@HEADER
// ************************************************************************
//
// containers_sparse_dense_matrix_product.hpp
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

#ifndef CONTAINERS_MATRIX_OPERATIONS_SPARSE_DENSE_MATRIX_PRODUCT_HPP_
#define CONTAINERS_MATRIX_OPERATIONS_SPARSE_DENSE_MATRIX_PRODUCT_HPP_

#include "../../meta/containers_matrix_meta.hpp"
#include "../concrete/containers_matrix_dense_sharedmem_eigen.hpp"
#include "../concrete/containers_matrix_sparse_sharedmem_eigen.hpp"
#include <EpetraExt_MatrixMatrix.h>
#include "TpetraExt_MatrixMatrix.hpp"
#include "Epetra_LocalMap.h"
#include <Epetra_Export.h>
#include <Epetra_Import.h>

namespace pressio{
namespace containers{
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
} // end namespace containers
}//end namespace pressio
#endif
