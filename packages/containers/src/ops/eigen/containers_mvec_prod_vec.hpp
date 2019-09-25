#ifndef CONTAINERS_SRC_OPS_EIGEN_MULTI_VECTOR_PROD_VECTOR_HPP_
#define CONTAINERS_SRC_OPS_EIGEN_MULTI_VECTOR_PROD_VECTOR_HPP_

#include "../containers_ops_meta.hpp"
#include "../../multi_vector/containers_multi_vector_meta.hpp"

namespace pressio{ namespace containers{ namespace ops{

/*
 * multi_vector prod vector
 */

template <
  typename mvec_type,
  typename vec_type,
  ::pressio::mpl::enable_if_t<
   containers::meta::is_multi_vector_wrapper_eigen<mvec_type>::value and
   containers::meta::is_vector_wrapper_eigen<vec_type>::value and
   containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value
    > * = nullptr
  >
void product(const mvec_type & mvA, const vec_type & vecB, vec_type & C){

  assert( C.size() == mvA.length() );
  //zero out result
  C.setZero();

  const auto numVecs = mvA.numVectors();
  // const auto Alength = mvA.length();
  // size of vecB
  assert(numVecs == vecB.size());

  // compute
  (*C.data()) = (*mvA.data()) * (*vecB.data());
  // for (size_t i=0; i<(size_t)Alength; i++){
  //   for (size_t j=0; j<(size_t)numVecs; j++){
  //     C[i] += mvA(i,j) * vecB[j];
  //   }
  // }
}//end function


// result is constructed and returned
template <
  typename mvec_type,
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_eigen<mvec_type>::value and
    containers::meta::is_vector_wrapper_eigen<vec_type>::value and
    containers::details::traits<vec_type>::is_dynamic and
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value
    > * = nullptr
  >
vec_type product(const mvec_type & mvA, const vec_type & vecB){

  vec_type c(mvA.length());
  product(mvA, vecB, c);
  return c;
}


}}}//end namespace pressio::containers::ops
#endif
