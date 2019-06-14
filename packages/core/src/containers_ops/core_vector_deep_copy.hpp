
#ifndef CORE_CONTAINER_OPS_VECTOR_DEEP_COPY_HPP_
#define CORE_CONTAINER_OPS_VECTOR_DEEP_COPY_HPP_

#include "core_ops_meta.hpp"
#include "../vector/core_vector_meta.hpp"
#ifdef HAVE_PYBIND11
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#endif

namespace rompp{ namespace core{ namespace ops{

//--------------------------------------------------------------------------
// enable for wrappers
//--------------------------------------------------------------------------
template<
  typename T,
  ::rompp::mpl::enable_if_t<
    ::rompp::core::meta::is_core_vector_wrapper<T>::value
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
    ::rompp::core::meta::is_cstyle_array_pybind11<T>::value
    > * = nullptr
  >
void deep_copy(const T & src, T & dest){
  if (src.ndim() > 1){
    throw std::runtime_error("core::ops::deep_copy: v.ndims()!=1. this operation currently supported for vectors only");
  }

  const auto vsz = src.size();
  if (vsz != dest.size())
    throw std::runtime_error("core::ops::deep_copy: Input shapes must match");

  for (decltype(dest.size()) i=0; i<vsz; ++i){
    dest.mutable_at(i) = src.at(i);
  }
}
#endif


}}}//end namespace rompp::core::ops
#endif
