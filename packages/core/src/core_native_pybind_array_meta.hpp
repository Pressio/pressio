
#ifdef HAVE_PYBIND11
#ifndef CORE_NATIVE_PYBIND11_ARRAY_HPP_
#define CORE_NATIVE_PYBIND11_ARRAY_HPP_

#include "./meta/core_meta_basic.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>


namespace rompp{ namespace core{ namespace meta {

/*
 * this metafunction is here because a pybind11::array_t
 * can have arbitrary size sine it maps to numpy.
 * And I don't know yet if a pybind11::array_t can be
 * checked at compile time to be a vector or matrix
 */

template <typename T, typename enable = void>
struct is_cstyle_array_pybind11 : std::false_type {};

template <typename T>
struct is_cstyle_array_pybind11<
  T,
  ::rompp::mpl::enable_if_t<
    mpl::is_same<
      T,
      pybind11::array_t<
	typename T::value_type,
	pybind11::array::c_style
	>
      >::value
    >
  > : std::true_type{};
//----------------------------------------------

template <typename T, typename enable = void>
struct is_array_pybind11 : std::false_type {};

template <typename T>
struct is_array_pybind11<
  T,
  ::rompp::mpl::enable_if_t<
    is_cstyle_array_pybind11<T>::value
    >
  > : std::true_type{};


}}}//end namespace rompp::core::meta
#endif
#endif
