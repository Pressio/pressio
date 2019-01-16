
#ifndef CORE_EIGEN_VECTOR_DOT_EIGEN_VECTOR_HPP_
#define CORE_EIGEN_VECTOR_DOT_EIGEN_VECTOR_HPP_

#include "../../core_ops_meta.hpp"
#include "../../../vector/core_vector_meta.hpp"

namespace rompp{ namespace core{ namespace ops{

// Eigen vector dot eigen vector
template <
  typename vec_type,
  core::meta::enable_if_t<
    core::meta::is_eigen_vector_wrapper<vec_type>::value
    > * = nullptr
  >
void dot(const vec_type & vecA,
	 const vec_type & vecB,
	 typename details::traits<vec_type>::scalar_t & result)
{
  assert(vecA.size() == vecB.size());
  result = vecA.data()->dot(*vecB.data());
}

// Eigen vector dot eigen vector, return value
template <
  typename vec_type,
  core::meta::enable_if_t<
    core::meta::is_eigen_vector_wrapper<vec_type>::value
    > * = nullptr
  >
auto dot(const vec_type & vecA,
	 const vec_type & vecB)
  -> typename details::traits<vec_type>::scalar_t
{
  using sc_t = typename details::traits<vec_type>::scalar_t;
  sc_t result ={};
  dot(vecA, vecB, result);
  return result;
}

}}} // end namespace rompp::core::ops
#endif
