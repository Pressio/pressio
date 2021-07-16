
#include <gtest/gtest.h>
#include "pressio_containers.hpp"

using scalar_t  = double;
using py_f_t	= pybind11::array_t<scalar_t, pybind11::array::f_style>;
using py_c_t	= pybind11::array_t<scalar_t, pybind11::array::c_style>;
using v_f_t	= pressio::containers::Vector<py_f_t>;
using v_c_t	= pressio::containers::Vector<py_c_t>;

TEST(containers_vector_pybind, native_meta){
  static_assert(!pressio::containers::predicates::is_cstyle_array_pybind<py_f_t>::value, "");
  static_assert(pressio::containers::predicates::is_fstyle_array_pybind<py_f_t>::value, "");
  static_assert(pressio::containers::predicates::is_array_pybind<py_f_t>::value, "");

  static_assert(pressio::containers::predicates::is_cstyle_array_pybind<py_c_t>::value, "");
  static_assert(!pressio::containers::predicates::is_fstyle_array_pybind<py_c_t>::value, "");
  static_assert(pressio::containers::predicates::is_array_pybind<py_c_t>::value, "");
}

TEST(containers_vector_pybind, vf_meta){
  static_assert(!pressio::containers::predicates::is_cstyle_array_pybind<v_f_t>::value, "");
  static_assert(!pressio::containers::predicates::is_fstyle_array_pybind<v_f_t>::value, "");
  static_assert(!pressio::containers::predicates::is_array_pybind<v_f_t>::value, "");
  static_assert(pressio::containers::predicates::is_fstyle_vector_wrapper_pybind<v_f_t>::value, "");
  static_assert(!pressio::containers::predicates::is_cstyle_vector_wrapper_pybind<v_f_t>::value, "");
  static_assert(pressio::containers::predicates::is_vector_wrapper_pybind<v_f_t>::value, "");
}

TEST(containers_vector_pybind, vc_meta){
  static_assert(!pressio::containers::predicates::is_cstyle_array_pybind<v_c_t>::value, "");
  static_assert(!pressio::containers::predicates::is_fstyle_array_pybind<v_c_t>::value, "");
  static_assert(!pressio::containers::predicates::is_array_pybind<v_c_t>::value, "");
  static_assert(!pressio::containers::predicates::is_fstyle_vector_wrapper_pybind<v_c_t>::value, "");
  static_assert(pressio::containers::predicates::is_cstyle_vector_wrapper_pybind<v_c_t>::value, "");
  static_assert(pressio::containers::predicates::is_vector_wrapper_pybind<v_c_t>::value, "");
}

TEST(containers_vector_pybind, f_traits){
  using traits = pressio::containers::details::traits<v_f_t>;
  static_assert(traits::wrapped_vector_identifier
		== pressio::containers::details::WrappedVectorIdentifier::Pybind, "");

  static_assert(std::is_same<typename traits::scalar_t, scalar_t>::value, "");
  static_assert(traits::is_static == false, "");
}
