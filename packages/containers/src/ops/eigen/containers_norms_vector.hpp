
#ifndef CONTAINERS_SRC_OPS_EIGEN_NORMS_HPP_
#define CONTAINERS_SRC_OPS_EIGEN_NORMS_HPP_

#include "../containers_ops_meta.hpp"
#include "../../vector/containers_vector_meta.hpp"

namespace pressio{ namespace containers{ namespace ops{

template <
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_eigen<vec_type>::value
    > * = nullptr
  >
auto norm1(const vec_type & a)
  -> typename details::traits<vec_type>::scalar_t
{
  using sc_t = typename details::traits<vec_type>::scalar_t;
  sc_t result = 0.0;
  for (decltype(a.size()) i=0; i<a.size(); i++)
    result += std::abs(a(i));
  return result;
}


template <
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_eigen<vec_type>::value
    > * = nullptr
  >
auto norm2(const vec_type & a)
  -> typename details::traits<vec_type>::scalar_t
{
  using sc_t = typename details::traits<vec_type>::scalar_t;
  sc_t result = 0.0;
  for (decltype(a.size()) i=0; i<a.size(); i++)
    result += a[i]*a[i];
  return std::sqrt(result);
}

}}}//end namespace pressio::containers::ops
#endif
