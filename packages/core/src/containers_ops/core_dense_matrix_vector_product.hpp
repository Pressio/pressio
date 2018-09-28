
#ifndef CORE_MATRIX_OPERATIONS_DENSE_MATRIX_VECTOR_PRODUCT_HPP_
#define CORE_MATRIX_OPERATIONS_DENSE_MATRIX_VECTOR_PRODUCT_HPP_

#include "core_ops_meta.hpp"
#include "../vector/core_vector_meta.hpp"
#include "../matrix/core_matrix_meta.hpp"

namespace rompp{
namespace core{
namespace ops{
  
/*---------------------------------------------------------
c = A b
- A is dense matrix from eigen
- b is vector from eigen
---------------------------------------------------------*/
  
template <typename mat_t, typename vec1_t, typename vec2_t,
	  core::meta::enable_if_t<
	    core::meta::is_eigen_dense_matrix_wrapper<mat_t>::value &&
	    core::meta::is_eigen_vector_wrapper<vec1_t>::value &&
	    core::meta::is_eigen_vector_wrapper<vec2_t>::value
	    > * = nullptr
	  >
void product(const mat_t & A, const vec1_t & b, vec2_t & c){

  assert(A.cols() == b.size());
  assert(c.size() == A.rows());
  (*c.data()) = (*A.data()) * (*b.data());
}

  
} // end namespace ops
} // end namespace core
}//end namespace rompp

#endif













// /*--------------------------------------------------------
//   EPETRA 
//   c = A b , 
//   - A = DENSE matrix 
//   - b = SINGLE vector
// -----------------------------------------------------------*/
  
// template <typename matrix_type,
// 	  typename vector_type,
// 	  core::meta::enable_if_t<
// 	    core::details::traits<matrix_type>::is_matrix==1 &&
// 	    core::details::traits<matrix_type>::isEpetra==1 &&
// 	    core::details::traits<matrix_type>::is_dense==1 &&
// 	    core::details::traits<vector_type>::is_vector==1 &&
// 	    core::details::traits<vector_type>::isEpetra==1
// 	    > * = nullptr>
// auto product(const matrix_type & A,
// 	     const vector_type & b)
// {

//    // I tried here to use the Multiply method of MultiVectors 
//    //   but it does not seem to work as expected. 
//    //   When A,b are all distributed, I don't get 
//    //   the right result. So we need to figure out why. 

//    //   Only solution that worked is to do this trick: 
//    //      b is distributed -> import into b replicated -> do multiply

//    //   Here we are doing matrix-vector product, where vector
//    //   is a single vector. So for now we replicate the
//    //   distributed vector across all processes.  
//    //   This is not too bad, but we should find out why not working 
//    //   for fully distributed case
    

//   const auto bGSize = b.globalSize();
//   assert( A.globalCols() == bGSize );
//   vector_type c( A.getDataMap() );
  
//   if ( b.isDistributedGlobally() ){
//     // define local map
//     Epetra_LocalMap locMap( bGSize, 0, b.commCRef() );
//     // define replicated vector
//     Epetra_Vector bRep(locMap);    
//     // get distributed map
//     auto & srcMap = b.getDataMap();
//     // define importer: Epetra_Import(targetMap, sourceMap)
//     Epetra_Import globToLocalImporter(locMap, srcMap);
//     // import global -> local
//     bRep.Import(*b.data(), globToLocalImporter, Insert);  
//     c.data()->Multiply( 'N','N', 1.0,  *A.data(), bRep, 0.0 );
//   }
//   else{
//     c.data()->Multiply( 'N','N', 1.0,  *A.data(), *b.data(), 0.0 );
//   }

//   return c;
// }

