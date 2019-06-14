
#ifndef CORE_CONTAINER_OPS_VECTOR_SET_ZERO_HPP_
#define CORE_CONTAINER_OPS_VECTOR_SET_ZERO_HPP_

#include "core_ops_meta.hpp"
#include "../vector/core_vector_meta.hpp"
#ifdef HAVE_PYBIND11
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#endif

//----------------------------------------------------------------------
//  overloads for setting all entries of vector to zero
//----------------------------------------------------------------------

namespace rompp{ namespace core{ namespace ops{

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
void set_zero(T & v){
  // make sure this is a vector
  if (v.ndim() > 1){
    throw std::runtime_error("core::ops::set_zero: v.ndims()!=1, while this operation requires a vector");
  }

  using scalar_t = typename T::value_type;
  for (decltype(v.size()) i=0; i<v.size(); ++i){
    v.mutable_at(i) = ::rompp::core::constants::zero<scalar_t>();
  }
}
#endif

}}}//end namespace rompp::core::ops
#endif
