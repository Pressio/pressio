
#ifndef CONTAINERS_CONTAINER_OPS_VECTOR_TWO_TERMS_UPDATE_HPP_
#define CONTAINERS_CONTAINER_OPS_VECTOR_TWO_TERMS_UPDATE_HPP_

#include "../containers_ops_meta.hpp"
#include "../../vector/containers_vector_meta.hpp"

#ifdef HAVE_TRILINOS
#include "containers_vector_do_update_kokkos_functors.hpp"
#endif

#ifdef HAVE_PYBIND11
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#endif

//----------------------------------------------------------------------
//  overloads for computing this: V = a * V + b * V1 + c * V2
//----------------------------------------------------------------------

namespace pressio{ namespace containers{ namespace ops{


//--------------------------------------------------------------------------
// enable for vectors supporting expression templates
//--------------------------------------------------------------------------
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper<T>::value and
    ::pressio::containers::meta::has_expression_templates_support<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t &a,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c){
  v = a*v + b*v1 + c*v2;
}

template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper<T>::value and
    ::pressio::containers::meta::has_expression_templates_support<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c){
  v = b*v1 + c*v2;
}


//--------------------------------------------------------------------------
// enable for pybind11::array_t
//--------------------------------------------------------------------------
#ifdef HAVE_PYBIND11
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_cstyle_array_pybind11<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t &a,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c){
  std::cout << " do_update pybind11 two terms " << std::endl;
  // make sure this is a vector
  if (v.ndim() > 1){
    throw std::runtime_error("containers::ops::do_update: v.ndims()!=1, while this operation requires a vector");
  }
  const auto vsz = v.size();
  if (vsz != v1.size() and vsz != v2.size())
    throw std::runtime_error("containers::ops::do_update: Input shapes must match");

  for (decltype(v.size()) i=0; i<vsz; ++i){
    v.mutable_at(i) = a*v.at(i) + b*v1.at(i) + c*v2.at(i);
  }
}

template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_cstyle_array_pybind11<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c){
  // make sure this is a vector
  if (v.ndim() > 1){
    throw std::runtime_error("containers::ops::do_update: v.ndims()!=1, while this operation requires a vector");
  }

  const auto vsz = v.size();
  if (vsz != v1.size() and vsz != v2.size())
    throw std::runtime_error("containers::ops::do_update: Input shapes must match");

  for (decltype(v.size()) i=0; i<vsz; ++i){
    v.mutable_at(i) = b*v1.at(i) + c*v2.at(i);
  }
}
#endif


//--------------------------------------------------------------------------
// enable for tpetra and tpetra block vectors
//--------------------------------------------------------------------------
#ifdef HAVE_TRILINOS
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_tpetra<T>::value or
    ::pressio::containers::meta::is_vector_wrapper_tpetra_block<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t &a,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c)
{
  constexpr auto one  = ::pressio::utils::constants::one<scalar_t>();
  v.data()->update(b, *v1.data(), a); // v = v + b * v1
  v.data()->update(c, *v2.data(), one); // add c*v2
}

template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_tpetra<T>::value or
    ::pressio::containers::meta::is_vector_wrapper_tpetra_block<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c)
{
  constexpr auto one  = ::pressio::utils::constants::one<scalar_t>();
  constexpr auto zero = ::pressio::utils::constants::zero<scalar_t>();

  v.data()->update(b, *v1.data(), zero); // v = b * v1
  v.data()->update(c, *v2.data(), one); // add c*v2
}
#endif


//--------------------------------------------------------------------------
// enable for Kokkos wrappers
//--------------------------------------------------------------------------
#ifdef HAVE_TRILINOS

template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_kokkos<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t &a,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c)
{
  using view_t = typename ::pressio::containers::details::traits<T>::wrapped_t;
  using fnctr_t = ::pressio::containers::ops::impl::DoUpdateTwoTermsFunctor<view_t, scalar_t>;
  fnctr_t F(*v.data(), *v1.data(), *v2.data(), a, b, c);
  Kokkos::parallel_for(v.size(), F);
}

template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_kokkos<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c)
{
  using view_t = typename ::pressio::containers::details::traits<T>::wrapped_t;
  using fnctr_t = ::pressio::containers::ops::impl::DoUpdateTwoTermsFunctor<view_t, scalar_t>;
  fnctr_t F(*v.data(), *v1.data(), *v2.data(), b, c);
  Kokkos::parallel_for(v.size(), F);
}

#endif

}}}//end namespace pressio::containers::ops
#endif
