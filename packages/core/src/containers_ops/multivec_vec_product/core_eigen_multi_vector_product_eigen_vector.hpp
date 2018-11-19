
#ifndef CORE_CONTAINER_OPS_MVEC_VEC_PROD_EIGEN_MULTI_VECTOR_PRODUCT_EIGEN_VECTOR_HPP_
#define CORE_CONTAINER_OPS_MVEC_VEC_PROD_EIGEN_MULTI_VECTOR_PRODUCT_EIGEN_VECTOR_HPP_

#include "../core_ops_meta.hpp"
#include "../../vector/core_vector_meta.hpp"
#include "../../multi_vector/core_multi_vector_meta.hpp"
#include "../../vector/concrete/core_vector_sharedmem_eigen_dynamic.hpp"
#ifdef HAVE_TRILINOS
#include "../../vector/concrete/core_vector_distributed_epetra.hpp"
#endif

namespace rompp{ namespace core{ namespace ops{

//-----------------------------------------------------
// Eigen multivector product with eigen vector
template <typename mvec_type,
	  typename vec_type,
  core::meta::enable_if_t<
   core::meta::is_eigen_multi_vector_wrapper<mvec_type>::value and
   core::meta::is_eigen_vector_wrapper<vec_type>::value and
   core::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value
    > * = nullptr
  >
void product(const mvec_type & mvA,
	     const vec_type & vecB,
	     core::Vector<Eigen::VectorXd> & C){

  assert( C.size() == mvA.length() );
  //zero out result
  C.setZero();
  // sizes of mvA
  auto numVecs = mvA.numVectors();
  auto Alength = mvA.length();

  // size of vecB
  size_t vecBLen = vecB.size();
  assert(size_t(numVecs) == vecBLen);

  // compute
  for (decltype(Alength) i=0; i<Alength; i++){
    for (decltype(numVecs) j=0; j<numVecs; j++){
      C[i] += mvA(i,j) * vecB[j];
    }
  }
}//end function


// result is constructed and returned
template <typename mvec_type,
	  typename vec_type,
  core::meta::enable_if_t<
   core::meta::is_eigen_multi_vector_wrapper<mvec_type>::value and
   core::meta::is_eigen_vector_wrapper<vec_type>::value and
   core::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value
    > * = nullptr
  >
auto product(const mvec_type & mvA, const vec_type & vecB)
-> core::Vector<
    Eigen::Matrix<typename core::details::traits<mvec_type>::scalar_t,
                  Eigen::Dynamic,1>
>{

  using sc_t = typename core::details::traits<mvec_type>::scalar_t;
  core::Vector<Eigen::Matrix<sc_t,Eigen::Dynamic,1>> c(mvA.length());
  product(mvA, vecB, c);
  return c;
}
//-------------------------------------------------------
//-------------------------------------------------------


}}}//end namespace rompp::core::ops
#endif
