
#ifndef CORE_MULTI_VECTOR_VECTOR_PRODUCT_HPP_
#define CORE_MULTI_VECTOR_VECTOR_PRODUCT_HPP_

#include "core_ops_meta.hpp"
#include "../vector/core_vector_meta.hpp"
#include "../multi_vector/core_multi_vector_meta.hpp"

namespace rompp{
namespace core{
namespace ops{


#ifdef HAVE_TRILINOS
//--------------------------------------------
//  Epetra multivector with eigen vector
//--------------------------------------------

// we pass the result object
template <typename mvec_type,
	  typename vec_type,
  core::meta::enable_if_t<
   core::meta::is_epetra_multi_vector_wrapper<mvec_type>::value and
   core::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value and
    (core::meta::is_eigen_vector_wrapper<vec_type>::value or
     core::meta::is_armadillo_column_vector_wrapper<vec_type>::value or
     core::meta::is_armadillo_row_vector_wrapper<vec_type>::value)
    > * = nullptr
  >
void product(const mvec_type & mvA,
	     const vec_type & vecB,
	     core::Vector<Epetra_Vector> & C){

  // using sc_t = typename details::traits<mvec_type>::scalar_t;
  // how many vectors are in mvA
  auto numVecs = mvA.globalNumVectors();
  // size of vecB
  auto vecBLen = vecB.size();
  assert(numVecs == vecBLen);
  
  // the data map of the multivector
  auto mvMap = mvA.getDataMap();
  // my number of rows
  auto myNrows = mvMap.NumMyElements();
  for (int i=0; i<myNrows; i++){
    for (decltype(numVecs) j=0; j<numVecs; j++){
      C[i] += mvA(i,j) * vecB[j];
    }
  }
}

// result is the return type
template <typename mvec_type,
	  typename vec_type,
  core::meta::enable_if_t<
   core::meta::is_epetra_multi_vector_wrapper<mvec_type>::value and
   core::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value and
    (core::meta::is_eigen_vector_wrapper<vec_type>::value or
     core::meta::is_armadillo_column_vector_wrapper<vec_type>::value or
     core::meta::is_armadillo_row_vector_wrapper<vec_type>::value)
  > * = nullptr
 >
auto product(const mvec_type & mvA,
	     const vec_type & vecB) {

  // here, mvA is distrubted, but vecB is NOT.
  // we interpret this as a linear combination of vectors
    
  // the data map of the multivector
  auto mvMap = mvA.getDataMap();
  // result is an Epetra Vector with same distribution of mvA  
  using res_t = core::Vector<Epetra_Vector>;
  res_t c(mvMap);
  // zero-out the result
  c.setZero();
  product(mvA, vecB, c);
  return c;
}
#endif //HAVE_TRILINOS



  
} // end namespace ops
} // end namespace core
}//end namespace rompp
#endif
