
#ifndef CORE_CONTAINER_OPS_VECTOR_FOUR_TERMS_UPDATE_HPP_
#define CORE_CONTAINER_OPS_VECTOR_FOUR_TERMS_UPDATE_HPP_

#include "../core_ops_meta.hpp"
#include "../../vector/core_vector_meta.hpp"

//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3 + e * V4
//----------------------------------------------------------------------

namespace rompp{ namespace core{ namespace ops{

// enable for vectors supporting expression templates
template<
  typename T,
  typename scalar_t,
  ::rompp::mpl::enable_if_t<
    ::rompp::core::meta::has_expression_templates_support<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t & a,
	       const T & v1, const scalar_t & b,
	       const T & v2, const scalar_t & c,
	       const T & v3, const scalar_t & d,
	       const T & v4, const scalar_t & e)
{
  v = a*v + b*v1 + c*v2 + d*v3 + e*v4;
}

template<
  typename T,
  typename scalar_t,
  ::rompp::mpl::enable_if_t<
    ::rompp::core::meta::has_expression_templates_support<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t & b,
	       const T & v2, const scalar_t & c,
	       const T & v3, const scalar_t & d,
	       const T & v4, const scalar_t & e)
{
  v = b*v1 + c*v2 + d*v3 + e*v4;
}


// enable for tpetra and tpetra block vectors NOT supporting expr templates
template<
  typename T,
  typename scalar_t,
  ::rompp::mpl::enable_if_t<
    ::rompp::core::meta::is_vector_wrapper_tpetra<T>::value or
    ::rompp::core::meta::is_vector_wrapper_tpetra_block<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t & a,
	       const T & v1, const scalar_t & b,
	       const T & v2, const scalar_t & c,
	       const T & v3, const scalar_t & d,
	       const T & v4, const scalar_t & e)
{
  constexpr auto one  = ::rompp::core::constants::one<scalar_t>();

  v.data()->update(b, *v1.data(), one); // v = v + b * v1
  v.data()->update(c, *v2.data(), one); // add c*v2
  v.data()->update(d, *v3.data(), one); // add d*v3
  v.data()->update(e, *v4.data(), one); // add e*v4
}

template<
  typename T,
  typename scalar_t,
  ::rompp::mpl::enable_if_t<
    ::rompp::core::meta::is_vector_wrapper_tpetra<T>::value or
    ::rompp::core::meta::is_vector_wrapper_tpetra_block<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t & b,
	       const T & v2, const scalar_t & c,
	       const T & v3, const scalar_t & d,
	       const T & v4, const scalar_t & e)
{
  constexpr auto one  = ::rompp::core::constants::one<scalar_t>();
  constexpr auto zero = ::rompp::core::constants::zero<scalar_t>();

  v.data()->update(b, *v1.data(), zero); // v = b * v1
  v.data()->update(c, *v2.data(), one); // add c*v2
  v.data()->update(d, *v3.data(), one); // add d*v3
  v.data()->update(e, *v4.data(), one); // add e*v4
}

}}}//end namespace rompp::core::ops
#endif
