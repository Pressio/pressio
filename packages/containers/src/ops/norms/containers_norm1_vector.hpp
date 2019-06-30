
#ifndef CONTAINERS_CONTAINER_OPS_NORMS_NORM1_VECTOR_HPP_
#define CONTAINERS_CONTAINER_OPS_NORMS_NORM1_VECTOR_HPP_

#include "../containers_ops_meta.hpp"
#include "../../vector/containers_vector_meta.hpp"

namespace rompp{ namespace containers{ namespace ops{

//--------------------------------------------------------
//  eigen vector wrapper
//--------------------------------------------------------
template <typename vec_type,
  ::rompp::mpl::enable_if_t<
    ::rompp::containers::meta::is_vector_wrapper_eigen<vec_type>::value
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

#ifdef HAVE_PYBIND11
// pybind11::array
template <typename vec_type,
  ::rompp::mpl::enable_if_t<
    ::rompp::containers::meta::is_array_pybind11<vec_type>::value
    > * = nullptr
  >
auto norm1(const vec_type & a)
  -> typename vec_type::value_type
{
  using sc_t = typename vec_type::value_type;
  sc_t result = 0.0;

  // make sure this is a vector
  if (a.ndim() > 1){
    throw std::runtime_error("a.ndims()!=1, this norm op is for pybind11 vectors");
  }

  for (decltype(a.size()) i=0; i<a.size(); i++)
    result += std::abs(a.at(i));
  return result;
}
#endif


//--------------------------------------------------------
//  teuchos serial dense vector wrapper
//--------------------------------------------------------
#ifdef HAVE_TRILINOS
template <typename vec_type,
  ::rompp::mpl::enable_if_t<
    ::rompp::containers::meta::is_dense_vector_wrapper_teuchos<vec_type>::value
    > * = nullptr
  >
auto norm1(const vec_type & a)
  -> typename details::traits<vec_type>::scalar_t
{
  using sc_t = typename details::traits<vec_type>::scalar_t;
  return static_cast<sc_t>(a.data()->normOne());
}
#endif 


}}}//end namespace rompp::containers::ops
#endif
