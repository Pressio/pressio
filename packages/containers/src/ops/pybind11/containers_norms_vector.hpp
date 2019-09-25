
#ifdef HAVE_PYBIND11
#ifndef CONTAINERS_SRC_OPS_PYBIND11_NORMS_HPP_
#define CONTAINERS_SRC_OPS_PYBIND11_NORMS_HPP_

#include "../containers_ops_meta.hpp"
#include "../../vector/containers_vector_meta.hpp"

namespace pressio{ namespace containers{ namespace ops{


template <
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_array_pybind11<vec_type>::value
    > * = nullptr
  >
auto norm1(const vec_type & a) -> typename vec_type::value_type
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


template <
  typename vec_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_array_pybind11<vec_type>::value
    > * = nullptr
  >
auto norm2(const vec_type & a) -> typename vec_type::value_type
{
  using sc_t = typename vec_type::value_type;
  sc_t result = 0.0;

  // make sure this is a vector
  if (a.ndim() > 1){
    throw std::runtime_error("a.ndims()!=1, this norm op is for pybind11 vectors");
  }

  for (decltype(a.size()) i=0; i<a.size(); i++)
    result += a.at(i)*a.at(i);
  return std::sqrt(result);
}


}}}//end namespace pressio::containers::ops
#endif
#endif
