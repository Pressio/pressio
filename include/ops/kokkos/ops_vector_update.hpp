/*
//@HEADER
// ************************************************************************
//
// ops_vector_update.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef OPS_KOKKOS_OPS_VECTOR_UPDATE_HPP_
#define OPS_KOKKOS_OPS_VECTOR_UPDATE_HPP_

#include <KokkosBlas1_axpby.hpp>

namespace pressio{ namespace ops{

namespace impl{
template <typename scalar_type, typename T1, typename ...Args>
struct _doUpdateAdmissibleOperands
{
  /* make sure we don't pass const objects to be modified.
     In kokkos it is legal to modify const views, not for pressio wrappers. */
  static_assert
    (!std::is_const<T1>::value,
     "ops:product: cannot modify a const-qualified wrapper of a Kokkos view");
  static_assert
    (containers::predicates::are_scalar_compatible<T1, Args...>::value and
     std::is_same<
     typename pressio::containers::details::traits<T1>::scalar_t,
     scalar_type>::value, "Types are not scalar compatible");
  static_assert
    (::pressio::containers::predicates::have_matching_execution_space<T1,Args...>::value,
     "operands need to have same execution space" );

  static constexpr auto value = true;
};
}//end namepsace ops::impl


//----------------------------------------------------------------------
// computing:  V = a * V + b * V1
//----------------------------------------------------------------------
template<typename T1, typename T2, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T1>::value and
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T2>::value
  >
update(T1 & v,	 const scalar_t & a,
	  const T2 & v1, const scalar_t & b)
{
  static_assert(impl::_doUpdateAdmissibleOperands<scalar_t, T1,T2>::value,"");
  ::KokkosBlas::axpby(b, *v1.data(), a, *v.data());
}

template<typename T1, typename T2, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T1>::value and
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T2>::value
  >
update(T1 & v, const T2 & v1, const scalar_t & b)
{
  static_assert(impl::_doUpdateAdmissibleOperands<scalar_t, T1,T2>::value,"");
  constexpr auto zero = ::pressio::utils::constants<scalar_t>::zero();
  ::KokkosBlas::axpby(b, *v1.data(), zero, *v.data());
}

//----------------------------------------------------------------------
//  overloads for computing this: V = a * V + b * V1 + c * V2
//----------------------------------------------------------------------
template<typename T1, typename T2, typename T3, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T1>::value and
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T2>::value and
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T3>::value
  >
update(T1 & v,	 const scalar_t &a,
	  const T2 & v1, const scalar_t &b,
	  const T3 & v2, const scalar_t &c)
{
  static_assert(impl::_doUpdateAdmissibleOperands<scalar_t, T1,T2,T3>::value,"");

  using view_t = typename ::pressio::containers::details::traits<T1>::wrapped_t;
  using fnctr_t = ::pressio::ops::impl::DoUpdateTwoTermsFunctor<view_t, scalar_t>;
  fnctr_t F(*v.data(), *v1.data(), *v2.data(), a, b, c);
  Kokkos::parallel_for(v.extent(0), F);
}

template<typename T1, typename T2, typename T3, typename scalar_t>
::pressio::mpl::enable_if_t<
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T1>::value and
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T2>::value and
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T3>::value
  >
update(T1 & v,
	  const T2 & v1, const scalar_t &b,
	  const T3 & v2, const scalar_t &c)
{
  static_assert(impl::_doUpdateAdmissibleOperands<scalar_t,T1,T2,T3>::value,"");

  using view_t = typename ::pressio::containers::details::traits<T1>::wrapped_t;
  using fnctr_t = ::pressio::ops::impl::DoUpdateTwoTermsFunctor<view_t, scalar_t>;
  fnctr_t F(*v.data(), *v1.data(), *v2.data(), b, c);
  Kokkos::parallel_for(v.extent(0), F);
}

//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3
//----------------------------------------------------------------------
template<
  typename T1, typename T2, typename T3, typename T4,
  typename scalar_t
  >
::pressio::mpl::enable_if_t<
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T1>::value and
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T2>::value and
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T3>::value and
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T4>::value
  >
update(T1 & v,	 const scalar_t &a,
	  const T2 & v1, const scalar_t &b,
	  const T3 & v2, const scalar_t &c,
	  const T4 & v3, const scalar_t &d)
{
  static_assert(impl::_doUpdateAdmissibleOperands<scalar_t, T1,T2,T3,T4>::value,"");

  using view_t = typename ::pressio::containers::details::traits<T1>::wrapped_t;
  using fnctr_t = ::pressio::ops::impl::DoUpdateThreeTermsFunctor<view_t, scalar_t>;
  fnctr_t F(*v.data(), *v1.data(), *v2.data(), *v3.data(), a, b, c, d);
  Kokkos::parallel_for(v.extent(0), F);
}

template<
  typename T1, typename T2, typename T3, typename T4,
  typename scalar_t
  >
::pressio::mpl::enable_if_t<
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T1>::value and
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T2>::value and
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T3>::value and
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T4>::value
  >
update(T1 & v,
	  const T2 & v1, const scalar_t &b,
	  const T3 & v2, const scalar_t &c,
	  const T4 & v3, const scalar_t &d)
{
  static_assert(impl::_doUpdateAdmissibleOperands<scalar_t, T1,T2,T3,T4>::value,"");

  using view_t = typename ::pressio::containers::details::traits<T1>::wrapped_t;
  using fnctr_t = ::pressio::ops::impl::DoUpdateThreeTermsFunctor<view_t, scalar_t>;
  fnctr_t F(*v.data(), *v1.data(), *v2.data(), *v3.data(), b, c, d);
  Kokkos::parallel_for(v.extent(0), F);
}

//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3 + e * V4
//----------------------------------------------------------------------
template<
  typename T1, typename T2, typename T3, typename T4, typename T5,
  typename scalar_t
  >
::pressio::mpl::enable_if_t<
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T1>::value and
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T2>::value and
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T3>::value and
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T4>::value
  >
update(T1 & v,	const scalar_t &a,
	  const T2 & v1, const scalar_t &b,
	  const T3 & v2, const scalar_t &c,
	  const T4 & v3, const scalar_t &d,
	  const T5 & v4, const scalar_t &e)
{
  static_assert(impl::_doUpdateAdmissibleOperands<scalar_t, T1,T2,T3,T4,T5>::value,"");

  using view_t = typename ::pressio::containers::details::traits<T1>::wrapped_t;
  using fnctr_t = ::pressio::ops::impl::DoUpdateFourTermsFunctor<view_t, scalar_t>;
  fnctr_t F(*v.data(), *v1.data(), *v2.data(), *v3.data(), *v4.data(), a, b, c, d, e);
  Kokkos::parallel_for(v.extent(0), F);
}

template<
  typename T1, typename T2, typename T3, typename T4, typename T5,
  typename scalar_t
  >
::pressio::mpl::enable_if_t<
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T1>::value and
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T2>::value and
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T3>::value and
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T4>::value
  >
update(T1 & v,
	  const T2 & v1, const scalar_t &b,
	  const T3 & v2, const scalar_t &c,
	  const T4 & v3, const scalar_t &d,
	  const T5 & v4, const scalar_t &e)
{
  static_assert(impl::_doUpdateAdmissibleOperands<scalar_t,T1,T2,T3,T4,T5>::value,"");

  using view_t = typename ::pressio::containers::details::traits<T1>::wrapped_t;
  using fnctr_t = ::pressio::ops::impl::DoUpdateFourTermsFunctor<view_t, scalar_t>;
  fnctr_t F(*v.data(), *v1.data(), *v2.data(), *v3.data(), *v4.data(), b, c, d, e);
  Kokkos::parallel_for(v.extent(0), F);
}

}}//end namespace pressio::ops
#endif  // OPS_KOKKOS_OPS_VECTOR_UPDATE_HPP_
