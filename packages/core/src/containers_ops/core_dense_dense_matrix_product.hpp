
#ifndef CORE_MATRIX_OPERATIONS_DENSE_DENSE_MATRIX_PRODUCT_HPP_
#define CORE_MATRIX_OPERATIONS_DENSE_DENSE_MATRIX_PRODUCT_HPP_

#include "core_ops_meta.hpp"
#include "../matrix/core_matrix_meta.hpp"

namespace rompp{
namespace core{
namespace ops{
  
/*---------------------------------------------------------
C = A B
- A is dense matrix from eigen (can be static or dynamic)
- B is dense matrix from eigen (can be static or dynamic)
---------------------------------------------------------*/
  
template <typename T1, typename T2, typename T3,
core::meta::enable_if_t<
  core::meta::is_eigen_dense_matrix_wrapper<T1>::value &&
  core::meta::is_eigen_dense_matrix_wrapper<T2>::value &&
  core::meta::is_eigen_dense_matrix_wrapper<T3>::value && 
  core::meta::wrapper_triplet_have_same_scalar<T1,T2,T3>::value 
  > * = nullptr
>
void product(const T1 & A, const T2 & B, T3 & C){

  assert(C.rows() == A.rows());
  assert(C.cols() == B.cols());
  (*C.data()) = (*A.data()) * (*B.data());
}

  
template <typename T1, typename T2,
     core::meta::enable_if_t<
       core::meta::is_eigen_dense_matrix_wrapper<T1>::value and
       core::meta::is_eigen_dense_matrix_wrapper<T2>::value and
       core::meta::wrapper_pair_have_same_scalar<T1,T2>::value
     > * = nullptr
    >
auto product(const T1 & A, const T2 & B){

  using sc_t = typename core::details::traits<T1>::scalar_t;
  using nat_t = Eigen::Matrix<sc_t,-1,-1>; 
  core::Matrix<nat_t> C(A.rows(),B.cols());
  product(A, B, C);
  return C;
}
  

  
} // end namespace ops
} // end namespace core
}//end namespace rompp
#endif





















// /*-----------------------------------------------------
//   C = A * B
//   A: epetra dense matrix
//   B: epetra dense matrix
//   C: epetra dense matrix
// -----------------------------------------------------*/
// template <typename mat_type,
//     typename std::enable_if<
//       details::traits<mat_type>::isEpetra &&
//       details::traits<mat_type>::is_dense
//       >::type * = nullptr
//     >
// auto product(const mat_type & A,
//        const mat_type & B)
// {

//   const auto BGSize = B.globalRows();
//   assert( A.globalCols() == BGSize );
  
//    // I tried here to use the Multiply method of MultiVectors 
//    //   but it does not seem to work as expected. 
//    //   When A,B are all distributed, I don't get 
//    //   the right result. So we need to figure out why. 

//    //   Only solution that worked is to do this trick: 
//    //      B is distributed -> import into B replicated -> do multiply

//    //   We should find out why not working for fully distributed case
    

//   // define local map
//   Epetra_LocalMap locMap( BGSize, 0, B.commCRef() );
//   // define replicated B
//   Epetra_MultiVector BRep(locMap, B.globalCols());
    
//   // get distributed map
//   auto & srcMap = B.getDataMap();
//   // define importer: Epetra_Import(targetMap, sourceMap)
//   Epetra_Import globToLocalImporter(locMap, srcMap);

//   // import global -> local
//   BRep.Import(*B.data(), globToLocalImporter, Insert);
  
//   mat_type C( A.getDataMap(), B.globalCols() );
//   C.data()->Multiply( 'N','N', 1.0,  *A.data(), BRep, 0.0 );  
//   return C;
// }
