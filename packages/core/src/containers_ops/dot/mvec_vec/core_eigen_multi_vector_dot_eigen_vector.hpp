
#ifndef CORE_EIGEN_MULTI_VECTOR_DOT_EIGEN_VECTOR_HPP_
#define CORE_EIGEN_MULTI_VECTOR_DOT_EIGEN_VECTOR_HPP_

#include "../../core_ops_meta.hpp"
#include "../../../vector/core_vector_meta.hpp"
#include "../../../multi_vector/core_multi_vector_meta.hpp"

namespace rompp{ namespace core{ namespace ops{

// Eigen multivector dot eigen vector
// result stored in Eigen DYNAMIC vector passed by reference
template <typename mvec_type,
	  typename vec_type,
	  typename result_vec_type,
  core::meta::enable_if_t<
    core::meta::is_eigen_multi_vector_wrapper<mvec_type>::value and
    core::meta::is_eigen_vector_wrapper<vec_type>::value and
    core::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value and
    core::meta::is_eigen_vector_wrapper<result_vec_type>::value and
    core::meta::wrapper_pair_have_same_scalar<vec_type,result_vec_type>::value and
    core::details::traits<result_vec_type>::is_dynamic
    > * = nullptr
  >
void dot(const mvec_type & mvA,
	 const vec_type & vecB,
	 result_vec_type & result){
  auto numVecs = mvA.numVectors();
  if ( result.size() != numVecs )
    result.resize(numVecs);
  *result.data() = (*mvA.data()).transpose() * (*vecB.data());
}
//--------------------------------------------------------


// Eigen multivector dot eigen vector
// result stored in Eigen STATIC vector passed by reference
template <typename mvec_type,
	  typename vec_type,
	  typename result_vec_type,
  core::meta::enable_if_t<
    core::meta::is_eigen_multi_vector_wrapper<mvec_type>::value and
    core::meta::is_eigen_vector_wrapper<vec_type>::value and
    core::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value and
    core::meta::is_eigen_vector_wrapper<result_vec_type>::value and
    core::meta::wrapper_pair_have_same_scalar<vec_type,result_vec_type>::value and
    core::details::traits<result_vec_type>::is_static
    > * = nullptr
  >
void dot(const mvec_type & mvA,
	 const vec_type & vecB,
	 result_vec_type & result){

  *result.data() = (*mvA.data()).transpose() * (*vecB.data());
}
//--------------------------------------------------------


// Eigen multivector dot eigen vector
// result is built and returned
template <typename mvec_type,
	  typename vec_type,
  core::meta::enable_if_t<
    core::meta::is_eigen_multi_vector_wrapper<mvec_type>::value and
    core::meta::is_eigen_vector_wrapper<vec_type>::value and
    core::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value
    > * = nullptr
  >
auto dot(const mvec_type & mvA, const vec_type & vecB)
-> core::Vector<
  Eigen::Matrix<typename core::details::traits<mvec_type>::scalar_t,
                Eigen::Dynamic, 1>
                >{

  using sc_t = typename core::details::traits<mvec_type>::scalar_t;
  core::Vector<Eigen::Matrix<sc_t, Eigen::Dynamic, 1>> c(vecB.size());
  dot(mvA,vecB,c);
  return c;
}
//--------------------------------------------------------

}}} // end namespace rompp::core::ops
#endif
