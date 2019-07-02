
#ifndef CONTAINERS_EIGEN_VECTOR_DOT_EIGEN_VECTOR_HPP_
#define CONTAINERS_EIGEN_VECTOR_DOT_EIGEN_VECTOR_HPP_

#include "../../containers_ops_meta.hpp"
#include "../../../vector/containers_vector_meta.hpp"

namespace pressio{ namespace containers{ namespace ops{

// Eigen vector dot eigen vector
template <
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_vector_wrapper_eigen<vec_type>::value
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
  ::pressio::mpl::enable_if_t<
    containers::meta::is_vector_wrapper_eigen<vec_type>::value
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

}}} // end namespace pressio::containers::ops
#endif
