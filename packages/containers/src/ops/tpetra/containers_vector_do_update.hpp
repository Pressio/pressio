
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_SRC_OPS_TPETRA_VECTOR_DO_UPDATE_HPP_
#define CONTAINERS_SRC_OPS_TPETRA_VECTOR_DO_UPDATE_HPP_

#include "../containers_ops_meta.hpp"
#include "../../vector/containers_vector_meta.hpp"

namespace pressio{ namespace containers{ namespace ops{

//----------------------------------------------------------------------
// computing:  V = a * V + b * V1
//----------------------------------------------------------------------
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_tpetra<T>::value
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
    ::pressio::containers::meta::is_vector_wrapper_tpetra<T>::value
    > * = nullptr
  >
void do_update(T & v, const T & v1, const scalar_t b)
{
  constexpr auto zero = ::pressio::utils::constants::zero<scalar_t>();
  v.data()->update(b, *v1.data(), zero); // v = b * v1
}


//----------------------------------------------------------------------
//  overloads for computing this: V = a * V + b * V1 + c * V2
//----------------------------------------------------------------------
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_tpetra<T>::value
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
    ::pressio::containers::meta::is_vector_wrapper_tpetra<T>::value
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


//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3
//----------------------------------------------------------------------
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_tpetra<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t &a,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d)
{
  constexpr auto one  = ::pressio::utils::constants::one<scalar_t>();

  v.data()->update(b, *v1.data(), a); // v = a*v + b*v1
  v.data()->update(c, *v2.data(), one); // add c*v2
  v.data()->update(d, *v3.data(), one); // add d*v3
}

template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_tpetra<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d)
{
  constexpr auto one  = ::pressio::utils::constants::one<scalar_t>();
  constexpr auto zero = ::pressio::utils::constants::zero<scalar_t>();

  v.data()->update(b, *v1.data(), zero); // v = b * v1
  v.data()->update(c, *v2.data(), one); // add c*v2
  v.data()->update(d, *v3.data(), one); // add d*v3
}


//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3 + e * V4
//----------------------------------------------------------------------
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_tpetra<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t &a,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d,
	       const T & v4, const scalar_t &e)
{
  constexpr auto one  = ::pressio::utils::constants::one<scalar_t>();

  v.data()->update(b, *v1.data(), a); // v = a*v + b*v1
  v.data()->update(c, *v2.data(), one); // add c*v2
  v.data()->update(d, *v3.data(), one); // add d*v3
  v.data()->update(e, *v4.data(), one); // add e*v4
}

template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_tpetra<T>::value
    > * = nullptr
  >
void do_update(T & v,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d,
	       const T & v4, const scalar_t &e)
{
  constexpr auto one  = ::pressio::utils::constants::one<scalar_t>();
  constexpr auto zero = ::pressio::utils::constants::zero<scalar_t>();

  v.data()->update(b, *v1.data(), zero); // v = b * v1
  v.data()->update(c, *v2.data(), one); // add c*v2
  v.data()->update(d, *v3.data(), one); // add d*v3
  v.data()->update(e, *v4.data(), one); // add e*v4
}


}}}//end namespace pressio::containers::ops
#endif
#endif
