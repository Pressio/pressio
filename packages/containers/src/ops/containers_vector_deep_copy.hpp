
#ifndef CONTAINERS_CONTAINER_OPS_VECTOR_DEEP_COPY_HPP_
#define CONTAINERS_CONTAINER_OPS_VECTOR_DEEP_COPY_HPP_

#include "containers_ops_meta.hpp"
#include "../vector/containers_vector_meta.hpp"
#ifdef HAVE_PYBIND11
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#endif

namespace rompp{ namespace containers{ namespace ops{

//--------------------------------------------------------------------------
// enable for wrappers
//--------------------------------------------------------------------------
template<
  typename T,
  ::rompp::mpl::enable_if_t<
    ::rompp::containers::meta::is_vector_wrapper<T>::value
    > * = nullptr
  >
void deep_copy(const T & src, T & dest){
  dest = src;
}

//--------------------------------------------------------------------------
// enable for pybind11::array_t
//--------------------------------------------------------------------------
#ifdef HAVE_PYBIND11
template<
  typename T,
  ::rompp::mpl::enable_if_t<
    ::rompp::containers::meta::is_cstyle_array_pybind11<T>::value
    > * = nullptr
  >
void deep_copy(const T & src, T & dest){
  if (src.ndim() > 1){
    throw std::runtime_error("containers::ops::deep_copy: v.ndims()!=1. this operation currently supported for vectors only");
  }

  const auto vsz = src.size();
  if (vsz != dest.size())
    throw std::runtime_error("containers::ops::deep_copy: Input shapes must match");

  for (decltype(dest.size()) i=0; i<vsz; ++i){
    dest.mutable_at(i) = src.at(i);
  }
}
#endif


}}}//end namespace rompp::containers::ops
#endif
