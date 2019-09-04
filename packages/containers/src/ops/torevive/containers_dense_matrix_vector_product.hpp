/*
//@HEADER
// ************************************************************************
//
// containers_dense_matrix_vector_product.hpp
//                     		      Pressio 
// Copyright 2019 National Technology & Engineering Solutions of Sandia,LLC 
//							      (NTESS)
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

#ifndef CONTAINERS_MATRIX_OPERATIONS_DENSE_MATRIX_VECTOR_PRODUCT_HPP_
#define CONTAINERS_MATRIX_OPERATIONS_DENSE_MATRIX_VECTOR_PRODUCT_HPP_

#include "../../meta/containers_vector_meta.hpp"
#include "../../vector/concrete/containers_vector_sharedmem_eigen.hpp"
#include "../../meta/containers_matrix_meta.hpp"
#include "../concrete/containers_matrix_dense_sharedmem_eigen.hpp"
#include "../concrete/containers_matrix_sparse_sharedmem_eigen.hpp"
#include "Epetra_LocalMap.h"
#include <Epetra_Export.h>
#include <Epetra_Import.h>

namespace pressio{
namespace containers{
namespace mat_ops{
  
/*---------------------------------------------------------
c = A b
- A is matrix from eigen
- b is vector from eigen
---------------------------------------------------------*/
  
template <typename matrix_type,
	  typename vector_t_1,
	  typename vector_t_2,
	  typename std::enable_if<
	    details::traits<vector_t_1>::isEigen &&
	    details::traits<vector_t_1>::is_vector &&
	    details::traits<vector_t_2>::isEigen &&
	    details::traits<vector_t_2>::is_vector &&
	    details::traits<matrix_type>::isEigen &&
	    details::traits<matrix_type>::is_dense &&
	    std::is_same<typename details::traits<matrix_type>::scalar_t,
			 typename details::traits<vector_t_1>::scalar_t
			 >::value &&
	    std::is_same<typename details::traits<matrix_type>::scalar_t,
			 typename details::traits<vector_t_2>::scalar_t
			 >::value 			 
	    >::type * = nullptr
	  >
void product(const matrix_type & A,
	     const vector_t_1 & b,
	     vector_t_2 & c)
{

  assert(A.cols() == b.size());
  assert(c.size() == A.rows());
  (*c.data()) = (*A.data()) * (*b.data());
}

  
  
/*--------------------------------------------------------
  EPETRA 
  c = A b , 
  - A = DENSE matrix 
  - b = SINGLE vector
-----------------------------------------------------------*/
  
template <typename matrix_type,
	  typename vector_type,
	  ::pressio::mpl::enable_if_t<
	    containers::details::traits<matrix_type>::is_matrix==1 &&
	    containers::details::traits<matrix_type>::isEpetra==1 &&
	    containers::details::traits<matrix_type>::is_dense==1 &&
	    containers::details::traits<vector_type>::is_vector==1 &&
	    containers::details::traits<vector_type>::isEpetra==1
	    > * = nullptr>
auto product(const matrix_type & A,
	     const vector_type & b)
{

  /* I tried here to use the Multiply method of MultiVectors 
     but it does not seem to work as expected. 
     When A,b are all distributed, I don't get 
     the right result. So we need to figure out why. 

     Only solution that worked is to do this trick: 
        b is distributed -> import into b replicated -> do multiply

     Here we are doing matrix-vector product, where vector
     is a single vector. So for now we replicate the
     distributed vector across all processes.  
     This is not too bad, but we should find out why not working 
     for fully distributed case
  */  

  const auto bGSize = b.globalSize();
  assert( A.globalCols() == bGSize );
  vector_type c( A.getDataMap() );
  
  if ( b.isDistributedGlobally() ){
    // define local map
    Epetra_LocalMap locMap( bGSize, 0, b.commCRef() );
    // define replicated vector
    Epetra_Vector bRep(locMap);    
    // get distributed map
    auto & srcMap = b.getDataMap();
    // define importer: Epetra_Import(targetMap, sourceMap)
    Epetra_Import globToLocalImporter(locMap, srcMap);
    // import global -> local
    bRep.Import(*b.data(), globToLocalImporter, Insert);  
    c.data()->Multiply( 'N','N', 1.0,  *A.data(), bRep, 0.0 );
  }
  else{
    c.data()->Multiply( 'N','N', 1.0,  *A.data(), *b.data(), 0.0 );
  }

  return c;
}
  

} // end namespace mat_ops
} // end namespace containers
}//end namespace pressio
#endif
