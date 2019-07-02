
#ifndef CONTAINERS_EIGEN_MULTI_VECTOR_DOT_EIGEN_VECTOR_HPP_
#define CONTAINERS_EIGEN_MULTI_VECTOR_DOT_EIGEN_VECTOR_HPP_

#include "../../containers_ops_meta.hpp"
#include "../../../vector/containers_vector_meta.hpp"
#include "../../../multi_vector/containers_multi_vector_meta.hpp"

namespace pressio{ namespace containers{ namespace ops{

// Eigen multivector dot eigen vector
// result stored in Eigen DYNAMIC vector passed by reference
template <typename mvec_type,
	  typename vec_type,
	  typename result_vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_eigen<mvec_type>::value and
    containers::meta::is_vector_wrapper_eigen<vec_type>::value and
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value and
    containers::meta::is_vector_wrapper_eigen<result_vec_type>::value and
    containers::meta::wrapper_pair_have_same_scalar<vec_type,result_vec_type>::value and
    containers::details::traits<result_vec_type>::is_dynamic
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
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_eigen<mvec_type>::value and
    containers::meta::is_vector_wrapper_eigen<vec_type>::value and
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value and
    containers::meta::is_vector_wrapper_eigen<result_vec_type>::value and
    containers::meta::wrapper_pair_have_same_scalar<vec_type,result_vec_type>::value and
    containers::details::traits<result_vec_type>::is_static
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
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_eigen<mvec_type>::value and
    containers::meta::is_vector_wrapper_eigen<vec_type>::value and
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value
    > * = nullptr
  >
auto dot(const mvec_type & mvA, const vec_type & vecB)
-> containers::Vector<
  Eigen::Matrix<typename containers::details::traits<mvec_type>::scalar_t,
                Eigen::Dynamic, 1>
                >{

  using sc_t = typename containers::details::traits<mvec_type>::scalar_t;
  containers::Vector<Eigen::Matrix<sc_t, Eigen::Dynamic, 1>> c(vecB.size());
  dot(mvA,vecB,c);
  return c;
}
//--------------------------------------------------------

}}} // end namespace pressio::containers::ops
#endif
