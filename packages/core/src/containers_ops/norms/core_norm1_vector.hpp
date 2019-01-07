
#ifndef CORE_CONTAINER_OPS_NORMS_NORM1_VECTOR_HPP_
#define CORE_CONTAINER_OPS_NORMS_NORM1_VECTOR_HPP_

#include "../core_ops_meta.hpp"
#include "../../vector/core_vector_meta.hpp"

namespace rompp{ namespace core{ namespace ops{

//--------------------------------------------------------
//  eigen vector wrapper
//--------------------------------------------------------
template <typename vec_type,
  core::meta::enable_if_t<
    core::meta::is_eigen_vector_wrapper<vec_type>::value
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
//--------------------------------------------------------

//--------------------------------------------------------
//  teuchos serial dense vector wrapper
//--------------------------------------------------------
#ifdef HAVE_TRILINOS
template <typename vec_type,
  core::meta::enable_if_t<
    core::meta::is_teuchos_serial_dense_vector_wrapper<vec_type>::value
    > * = nullptr
  >
auto norm1(const vec_type & a)
  -> typename details::traits<vec_type>::scalar_t
{
  using sc_t = typename details::traits<vec_type>::scalar_t;
  return static_cast<sc_t>(a.data()->normOne());
}
#endif 


}}}//end namespace rompp::core::ops
#endif
