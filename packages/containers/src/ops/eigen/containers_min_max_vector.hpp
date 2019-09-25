
#ifndef CONTAINERS_SRC_OPS_EIGEN_MINMAX_VECTOR_HPP_
#define CONTAINERS_SRC_OPS_EIGEN_MINMAX_VECTOR_HPP_

#include "../containers_ops_meta.hpp"
#include "../../vector/containers_vector_meta.hpp"

namespace pressio{ namespace containers{ namespace ops{


template <
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_eigen<vec_type>::value
    > * = nullptr
  >
auto max(const vec_type & a) -> typename details::traits<vec_type>::scalar_t
{
  return a.data()->maxCoeff();
}

template <
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_eigen<vec_type>::value
    > * = nullptr
  >
auto min(const vec_type & a) -> typename details::traits<vec_type>::scalar_t
{
  return a.data()->minCoeff();
}

}}}//end namespace pressio::containers::ops
#endif
