
#ifndef CONTAINERS_CONTAINER_OPS_VECTOR_ONE_TERM_UPDATE_HPP_
#define CONTAINERS_CONTAINER_OPS_VECTOR_ONE_TERM_UPDATE_HPP_

#include "../containers_ops_meta.hpp"
#include "../../vector/containers_vector_meta.hpp"
#ifdef HAVE_PYBIND11
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#endif
#ifdef HAVE_TRILINOS
#include<KokkosBlas1_axpby.hpp>
#endif

//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1
//----------------------------------------------------------------------

namespace pressio{ namespace containers{ namespace ops{


//----------------------------------------------------------------------
// enable for vectors supporting expression templates
//----------------------------------------------------------------------
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper<T>::value and
    ::pressio::containers::meta::has_expression_templates_support<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t a,
	       const T & v1, const scalar_t b)
{
  v = a*v + b*v1;
}

template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper<T>::value and
    ::pressio::containers::meta::has_expression_templates_support<T>::value
    > * = nullptr
  >
void do_update(T & v, const T & v1, const scalar_t  b)
{
  v = b*v1;
}


//--------------------------------------------------------------------
// enable for pybind11::array_t
//--------------------------------------------------------------------
#ifdef HAVE_PYBIND11
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_cstyle_array_pybind11<T>::value
    > * = nullptr
  >
void do_update(T & v, scalar_t a,
	       const T & v1, scalar_t b){
  // make sure this is a vector
  if (v.ndim() > 1){
    throw std::runtime_error("containers::ops::do_update: v.ndims()!=1, while this operation requires a vector");
  }

  const auto vsz = v.size();
  if (vsz != v1.size())
    throw std::runtime_error("containers::ops::do_update: Input shapes must match");

  for (decltype(v.size()) i=0; i<vsz; ++i){
    v.mutable_at(i) = a*v.at(i) + b*v1.at(i);
  }
}

template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_cstyle_array_pybind11<T>::value
    > * = nullptr
  >
void do_update(T & v, const T & v1, const scalar_t b){
  // make sure this is a vector
  if (v.ndim() > 1){
    throw std::runtime_error("containers::ops::do_update: v.ndims()!=1, while this operation requires a vector");
  }

  const auto vsz = v.size();
  if (vsz != v1.size())
    throw std::runtime_error("containers::ops::do_update: Input shapes must match");

  for (decltype(v.size()) i=0; i<vsz; ++i){
    v.mutable_at(i) = b*v1.at(i);
  }
}
#endif


//---------------------------------------------------------------------
// enable for tpetra and tpetra block vectors NOT supporting expr templates
//---------------------------------------------------------------------
#ifdef HAVE_TRILINOS
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_tpetra<T>::value or
    ::pressio::containers::meta::is_vector_wrapper_tpetra_block<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t a,
	       const T & v1, const scalar_t b)
{
  v.data()->update(b, *v1.data(), a); // v = a*v + b * v1
}

template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_tpetra<T>::value or
    ::pressio::containers::meta::is_vector_wrapper_tpetra_block<T>::value
    > * = nullptr
  >
void do_update(T & v, const T & v1, const scalar_t b)
{
  constexpr auto zero = ::pressio::utils::constants::zero<scalar_t>();
  v.data()->update(b, *v1.data(), zero); // v = b * v1
}
#endif


//---------------------------------------------------------------------
// enable for kokkos wrapper
//---------------------------------------------------------------------
#ifdef HAVE_TRILINOS
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_kokkos<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t & a,
	       const T & v1, const scalar_t & b)
{
  // v = a*v + b * v1
  KokkosBlas::axpby(b, *v1.data(), a, *v.data());
}

template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_kokkos<T>::value
    > * = nullptr
  >
void do_update(T & v, const T & v1, const scalar_t & b)
{
  // v = b*v1
  constexpr auto zero = ::pressio::utils::constants::zero<scalar_t>();
  KokkosBlas::axpby(b, *v1.data(), zero, *v.data());
}
#endif


}}}//end namespace pressio::containers::ops
#endif
