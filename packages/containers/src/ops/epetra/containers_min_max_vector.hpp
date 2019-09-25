
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_SRC_OPS_EPETRA_MINMAX_VECTOR_HPP_
#define CONTAINERS_SRC_OPS_EPETRA_MINMAX_VECTOR_HPP_

#include "../containers_ops_meta.hpp"
#include "../../vector/containers_vector_meta.hpp"

namespace pressio{ namespace containers{ namespace ops{


template <
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_epetra<vec_type>::value
    > * = nullptr
  >
auto max(const vec_type & a) -> typename details::traits<vec_type>::scalar_t
{
  typename details::traits<vec_type>::scalar_t result = {};
  a.data()->maxValue(&result);
  return result;
}

template <
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_epetra<vec_type>::value
    > * = nullptr
  >
auto min(const vec_type & a) -> typename details::traits<vec_type>::scalar_t
{
  typename details::traits<vec_type>::scalar_t result = {};
  a.data()->minValue(&result);
  return result;
}

}}}//end namespace pressio::containers::ops
#endif
#endif
