
#ifndef CONTAINERS_CONTAINER_OPS_MVEC_VEC_PROD_EIGEN_MULTI_VECTOR_PRODUCT_EIGEN_VECTOR_HPP_
#define CONTAINERS_CONTAINER_OPS_MVEC_VEC_PROD_EIGEN_MULTI_VECTOR_PRODUCT_EIGEN_VECTOR_HPP_

#include "../containers_ops_meta.hpp"
#include "../../vector/containers_vector_meta.hpp"
#include "../../multi_vector/containers_multi_vector_meta.hpp"
#include "../../vector/concrete/containers_vector_sharedmem_eigen_dynamic.hpp"
#ifdef HAVE_TRILINOS
#include "../../vector/concrete/containers_vector_distributed_epetra.hpp"
#endif

namespace pressio{ namespace containers{ namespace ops{

//-----------------------------------------------------
// Eigen multivector product with eigen vector
template <typename mvec_type,
	  typename vec_type,
  ::pressio::mpl::enable_if_t<
   containers::meta::is_multi_vector_wrapper_eigen<mvec_type>::value and
   containers::meta::is_vector_wrapper_eigen<vec_type>::value and
   containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value
    > * = nullptr
  >
void product(const mvec_type & mvA,
	     const vec_type & vecB,
	     containers::Vector<Eigen::VectorXd> & C){

  assert( C.size() == mvA.length() );
  //zero out result
  C.setZero();
  // sizes of mvA
  auto numVecs = mvA.numVectors();
  auto Alength = mvA.length();

  // size of vecB
  assert(size_t(numVecs) == size_t(vecB.size()));

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
  ::pressio::mpl::enable_if_t<
   containers::meta::is_multi_vector_wrapper_eigen<mvec_type>::value and
   containers::meta::is_vector_wrapper_eigen<vec_type>::value and
   containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value
    > * = nullptr
  >
auto product(const mvec_type & mvA, const vec_type & vecB)
-> containers::Vector<
    Eigen::Matrix<typename containers::details::traits<mvec_type>::scalar_t,
                  Eigen::Dynamic,1>
>{

  using sc_t = typename containers::details::traits<mvec_type>::scalar_t;
  containers::Vector<Eigen::Matrix<sc_t,Eigen::Dynamic,1>> c(mvA.length());
  product(mvA, vecB, c);
  return c;
}
//-------------------------------------------------------
//-------------------------------------------------------


}}}//end namespace pressio::containers::ops
#endif
