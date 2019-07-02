
#ifndef CONTAINERS_EIGEN_VECTOR_DOT_EIGEN_MULTI_VECTOR_HPP_
#define CONTAINERS_EIGEN_VECTOR_DOT_EIGEN_MULTI_VECTOR_HPP_

#include "../../containers_ops_meta.hpp"
#include "../../../vector/containers_vector_meta.hpp"
#include "../../../multi_vector/containers_multi_vector_meta.hpp"

namespace pressio{ namespace containers{ namespace ops{

// Eigen vector dot multivector
// result stored in Eigen DYNAMIC vector passed by reference
template <
  typename vec_type,
  typename mvec_type,
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
void dot(const vec_type & vec,
	 const mvec_type & mv,
	 result_vec_type & result)
{
  dot(mv, vec, result);
}


// Eigen vector dot multivector
// result stored in Eigen STATIC vector passed by reference
template <
  typename vec_type,
  typename mvec_type,
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
void dot(const vec_type & vec,
	 const mvec_type & mv,
	 result_vec_type & result){
  dot(mv, vec, result);
}


// Eigen vector dot multivector
// result is built and returned
template <typename vec_type,
	  typename mvec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_eigen<mvec_type>::value and
    containers::meta::is_vector_wrapper_eigen<vec_type>::value and
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value
    > * = nullptr
  >
auto dot(const vec_type & vec, const mvec_type & mv)
-> containers::Vector<
  Eigen::Matrix<typename containers::details::traits<mvec_type>::scalar_t,
                Eigen::Dynamic, 1>
                >{

  return dot(mv, vec);
}

}}} // end namespace pressio::containers::ops
#endif
