/*
//@HEADER
// ************************************************************************
//
// containers_dense_dense_matrix_product.hpp
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

#ifndef CONTAINERS_MATRIX_OPERATIONS_DENSE_DENSE_MATRIX_PRODUCT_HPP_
#define CONTAINERS_MATRIX_OPERATIONS_DENSE_DENSE_MATRIX_PRODUCT_HPP_

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
  
/*-----------------------------------------------------
  C = A * B

  A: epetra dense matrix
  B: epetra dense matrix
  C: epetra dense matrix
-----------------------------------------------------*/
template <typename mat_type,
	  typename std::enable_if<
	    details::traits<mat_type>::isEpetra &&
	    details::traits<mat_type>::is_dense
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
} // end namespace containers
}//end namespace pressio
#endif
