
#ifdef HAVE_KOKKOS
#ifndef CONTAINERS_SRC_OPS_KOKKOS_VECTOR_DO_UPDATE_HPP_
#define CONTAINERS_SRC_OPS_KOKKOS_VECTOR_DO_UPDATE_HPP_

#include "../containers_ops_meta.hpp"
#include "../../vector/containers_vector_meta.hpp"
#include "containers_vector_do_update_kokkos_functors.hpp"
#include<KokkosBlas1_axpby.hpp>

namespace pressio{ namespace containers{ namespace ops{

//----------------------------------------------------------------------
// computing:  V = a * V + b * V1
//----------------------------------------------------------------------
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


//----------------------------------------------------------------------
//  overloads for computing this: V = a * V + b * V1 + c * V2
//----------------------------------------------------------------------
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


//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3
//----------------------------------------------------------------------
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_kokkos<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t &a,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d)
{
  using view_t = typename ::pressio::containers::details::traits<T>::wrapped_t;
  using fnctr_t = ::pressio::containers::ops::impl::DoUpdateThreeTermsFunctor<view_t, scalar_t>;
  fnctr_t F(*v.data(), *v1.data(), *v2.data(), *v3.data(), a, b, c, d);
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
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d)
{
  using view_t = typename ::pressio::containers::details::traits<T>::wrapped_t;
  using fnctr_t = ::pressio::containers::ops::impl::DoUpdateThreeTermsFunctor<view_t, scalar_t>;
  fnctr_t F(*v.data(), *v1.data(), *v2.data(), *v3.data(), b, c, d);
  Kokkos::parallel_for(v.size(), F);
}


//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3 + e * V4
//----------------------------------------------------------------------
template<
  typename T,
  typename scalar_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_kokkos<T>::value
    > * = nullptr
  >
void do_update(T & v, const scalar_t &a,
	       const T & v1, const scalar_t &b,
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d,
	       const T & v4, const scalar_t &e)
{
  using view_t = typename ::pressio::containers::details::traits<T>::wrapped_t;
  using fnctr_t = ::pressio::containers::ops::impl::DoUpdateFourTermsFunctor<view_t, scalar_t>;
  fnctr_t F(*v.data(), *v1.data(), *v2.data(), *v3.data(), *v4.data(), a, b, c, d, e);
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
	       const T & v2, const scalar_t &c,
	       const T & v3, const scalar_t &d,
	       const T & v4, const scalar_t &e)
{
  using view_t = typename ::pressio::containers::details::traits<T>::wrapped_t;
  using fnctr_t = ::pressio::containers::ops::impl::DoUpdateFourTermsFunctor<view_t, scalar_t>;
  fnctr_t F(*v.data(), *v1.data(), *v2.data(), *v3.data(), *v4.data(), b, c, d, e);
  Kokkos::parallel_for(v.size(), F);
}


}}}//end namespace pressio::containers::ops
#endif
#endif
