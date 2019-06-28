
#ifndef ALGEBRA_CONTAINER_OPS_VECTOR_FOUR_TERMS_UPDATE_HPP_
#define ALGEBRA_CONTAINER_OPS_VECTOR_FOUR_TERMS_UPDATE_HPP_

#include "../algebra_ops_meta.hpp"
#include "../../vector/algebra_vector_meta.hpp"
#ifdef HAVE_PYBIND11
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#endif

//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3 + e * V4
//----------------------------------------------------------------------

namespace rompp{ namespace algebra{ namespace ops{


//---------------------------------------------------------------
// enable for vectors supporting expression templates
//---------------------------------------------------------------
template<
  typename T,
  typename scalar_t,
  ::rompp::mpl::enable_if_t<
    ::rompp::algebra::meta::is_algebra_vector_wrapper<T>::value and
    ::rompp::algebra::meta::has_expression_templates_support<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t a,
	       const T & v1, const scalar_t b,
	       const T & v2, const scalar_t c,
	       const T & v3, const scalar_t d,
	       const T & v4, const scalar_t e){
  v = a*v + b*v1 + c*v2 + d*v3 + e*v4;
}

template<
  typename T,
  typename scalar_t,
  ::rompp::mpl::enable_if_t<
    ::rompp::algebra::meta::is_algebra_vector_wrapper<T>::value and
    ::rompp::algebra::meta::has_expression_templates_support<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t b,
	       const T & v2, const scalar_t c,
	       const T & v3, const scalar_t d,
	       const T & v4, const scalar_t e){
  v = b*v1 + c*v2 + d*v3 + e*v4;
}



//--------------------------------------------------------------------------
// enable for pybind11::array_t
//--------------------------------------------------------------------------
#ifdef HAVE_PYBIND11
template<
  typename T,
  typename scalar_t,
  ::rompp::mpl::enable_if_t<
    ::rompp::algebra::meta::is_cstyle_array_pybind11<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t a,
	       const T & v1, const scalar_t b,
	       const T & v2, const scalar_t c,
	       const T & v3, const scalar_t d,
	       const T & v4, const scalar_t e){
  // make sure this is a vector
  if (v.ndim() > 1){
    throw std::runtime_error("algebra::ops::do_update: v.ndims()!=1, while this operation requires a vector");
  }

  const auto vsz = v.size();
  if (vsz != v1.size() and vsz != v2.size()
      and vsz != v3.size() and vsz != v4.size())
    throw std::runtime_error("algebra::ops::do_update: Input shapes must match");

  for (decltype(v.size()) i=0; i<vsz; ++i){
    v.mutable_at(i) = a*v.at(i) + b*v1.at(i) + c*v2.at(i) + d*v3.at(i) + e*v4.at(i);
  }
}

template<
  typename T,
  typename scalar_t,
  ::rompp::mpl::enable_if_t<
    ::rompp::algebra::meta::is_cstyle_array_pybind11<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t b,
	       const T & v2, const scalar_t c,
	       const T & v3, const scalar_t d,
	       const T & v4, const scalar_t e){
  // make sure this is a vector
  if (v.ndim() > 1){
    throw std::runtime_error("algebra::ops::do_update: v.ndims()!=1, while this operation requires a vector");
  }

  const auto vsz = v.size();
  if (vsz != v1.size() and vsz != v2.size()
      and vsz != v3.size() and vsz != v4.size())
    throw std::runtime_error("algebra::ops::do_update: Input shapes must match");

  for (decltype(v.size()) i=0; i<vsz; ++i){
    v.mutable_at(i) = b*v1.at(i) + c*v2.at(i) + d*v3.at(i) + e*v4.at(i);
  }
}
#endif



//-----------------------------------------------------------------------------
// enable for tpetra and tpetra block vectors
//-----------------------------------------------------------------------------
#ifdef HAVE_TRILINOS
template<
  typename T,
  typename scalar_t,
  ::rompp::mpl::enable_if_t<
    ::rompp::algebra::meta::is_vector_wrapper_tpetra<T>::value or
    ::rompp::algebra::meta::is_vector_wrapper_tpetra_block<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t a,
	       const T & v1, const scalar_t b,
	       const T & v2, const scalar_t c,
	       const T & v3, const scalar_t d,
	       const T & v4, const scalar_t e)
{
  constexpr auto one  = ::rompp::utils::constants::one<scalar_t>();

  v.data()->update(b, *v1.data(), a); // v = a*v + b*v1
  v.data()->update(c, *v2.data(), one); // add c*v2
  v.data()->update(d, *v3.data(), one); // add d*v3
  v.data()->update(e, *v4.data(), one); // add e*v4
}

template<
  typename T,
  typename scalar_t,
  ::rompp::mpl::enable_if_t<
    ::rompp::algebra::meta::is_vector_wrapper_tpetra<T>::value or
    ::rompp::algebra::meta::is_vector_wrapper_tpetra_block<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t b,
	       const T & v2, const scalar_t c,
	       const T & v3, const scalar_t d,
	       const T & v4, const scalar_t e)
{
  constexpr auto one  = ::rompp::utils::constants::one<scalar_t>();
  constexpr auto zero = ::rompp::utils::constants::zero<scalar_t>();

  v.data()->update(b, *v1.data(), zero); // v = b * v1
  v.data()->update(c, *v2.data(), one); // add c*v2
  v.data()->update(d, *v3.data(), one); // add d*v3
  v.data()->update(e, *v4.data(), one); // add e*v4
}
#endif

}}}//end namespace rompp::algebra::ops
#endif
