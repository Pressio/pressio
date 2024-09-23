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
#include "ops_vector_update_kokkos_functors.hpp"

namespace pressio{ namespace ops{

namespace impl{
template <typename scalar_type, typename T1, typename ...Args>
struct _kokkosUpdateAdmissibleOperands
{
  /* make sure we don't pass const objects to be modified.
     In kokkos it is legal to modify const views, not for pressio wrappers. */
  static_assert
    (!std::is_const<T1>::value,
     "ops:update: cannot modify a const-qualified wrapper of a Kokkos view");

  static constexpr auto value = true;
};

}//end namepsace ops::impl

//----------------------------------------------------------------------
// computing:  V = a * V + b * V1
//----------------------------------------------------------------------
template<
  typename T, typename T1,
  typename a_Type, typename b_Type
  >
std::enable_if_t<
  // rank-1 update common constraints
     ::pressio::Traits<T>::rank == 1
  && ::pressio::Traits<T1>::rank == 1
  // TPL/container specific
  && (::pressio::is_native_container_kokkos<T>::value
   || ::pressio::is_expression_acting_on_kokkos<T>::value)
  && (::pressio::is_native_container_kokkos<T1>::value
   ||  ::pressio::is_expression_acting_on_kokkos<T1>::value)
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<T, T1>::value
  && (std::is_floating_point<typename ::pressio::Traits<T>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<T>::scalar_type>::value)
  && std::is_convertible<a_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<b_Type, typename ::pressio::Traits<T>::scalar_type>::value
  >
update(T & v, const a_Type &a,
    const T1 & v1, const b_Type &b)
{
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v1, 0));

  using scalar_t = typename ::pressio::Traits<T>::scalar_type;
  static_assert(impl::_kokkosUpdateAdmissibleOperands<scalar_t,T,T1>::value,"");
  const scalar_t a_(a);
  const scalar_t b_(b);
  ::KokkosBlas::axpby(b_, impl::get_native(v1), a_, impl::get_native(v));
}

template<typename T, typename T1, typename b_Type>
std::enable_if_t<
  // rank-1 update common constraints
     ::pressio::Traits<T>::rank == 1
  && ::pressio::Traits<T1>::rank == 1
  // TPL/container specific
  && (::pressio::is_native_container_kokkos<T>::value
   || ::pressio::is_expression_acting_on_kokkos<T>::value)
  && (::pressio::is_native_container_kokkos<T1>::value
   ||  ::pressio::is_expression_acting_on_kokkos<T1>::value)
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<T, T1>::value
  && (std::is_floating_point<typename ::pressio::Traits<T>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<T>::scalar_type>::value)
  && std::is_convertible<b_Type, typename ::pressio::Traits<T>::scalar_type>::value
  >
update(T & v, const T1 & v1, const b_Type &b)
{
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v1, 0));

  using scalar_t = typename ::pressio::Traits<T>::scalar_type;
  static_assert(impl::_kokkosUpdateAdmissibleOperands<scalar_t,T,T1>::value,"");
  const scalar_t b_(b);
  ::KokkosBlas::axpby(b_, impl::get_native(v1), static_cast<scalar_t>(0), impl::get_native(v));
}

//----------------------------------------------------------------------
//  overloads for computing this: V = a * V + b * V1 + c * V2
//----------------------------------------------------------------------
template<
  typename T, typename T1, typename T2,
  typename a_Type, typename b_Type, typename c_Type
  >
std::enable_if_t<
  // rank-1 update common constraints
     ::pressio::Traits<T>::rank == 1
  && ::pressio::Traits<T1>::rank == 1
  && ::pressio::Traits<T2>::rank == 1
  // TPL/container specific
  && (::pressio::is_native_container_kokkos<T>::value
   || ::pressio::is_expression_acting_on_kokkos<T>::value)
  && (::pressio::is_native_container_kokkos<T1>::value
   ||  ::pressio::is_expression_acting_on_kokkos<T1>::value)
  && (::pressio::is_native_container_kokkos<T2>::value
   ||  ::pressio::is_expression_acting_on_kokkos<T2>::value)
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<T, T1, T2>::value
  && (std::is_floating_point<typename ::pressio::Traits<T>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<T>::scalar_type>::value)
  && std::is_convertible<a_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<b_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<c_Type, typename ::pressio::Traits<T>::scalar_type>::value
  >
update(T & v, const a_Type &a,
    const T1 & v1, const b_Type &b,
    const T2 & v2, const c_Type &c)
{
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v1, 0));
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v2, 0));

  using scalar_t = typename ::pressio::Traits<T>::scalar_type;
  static_assert(impl::_kokkosUpdateAdmissibleOperands<scalar_t,T,T1,T2>::value,"");

  const scalar_t a_(a);
  const scalar_t b_(b);
  const scalar_t c_(c);

  constexpr auto zero = static_cast<scalar_t>(0);
  if (b_ == zero) {
    ::pressio::ops::update(v, a_, v2, c_);
  } else if (c_ == zero) {
    ::pressio::ops::update(v, a_, v1, b_);
  } else {
    using v_t = typename impl::NativeType<T>::type;
    using v1_t = typename impl::NativeType<T1>::type;
    using v2_t = typename impl::NativeType<T2>::type;

    using fnctr_t = ::pressio::ops::impl::DoUpdateTwoTermsFunctor<v_t,v1_t,v2_t,scalar_t>;
    fnctr_t F(impl::get_native(v),
              impl::get_native(v1),
              impl::get_native(v2), a_, b_, c_);
    Kokkos::parallel_for(v.extent(0), F);
  }
}

template<
  typename T, typename T1, typename T2,
  typename b_Type, typename c_Type
  >
std::enable_if_t<
  // rank-1 update common constraints
     ::pressio::Traits<T>::rank == 1
  && ::pressio::Traits<T1>::rank == 1
  && ::pressio::Traits<T2>::rank == 1
  // TPL/container specific
  && (::pressio::is_native_container_kokkos<T>::value
   || ::pressio::is_expression_acting_on_kokkos<T>::value)
  && (::pressio::is_native_container_kokkos<T1>::value
   ||  ::pressio::is_expression_acting_on_kokkos<T1>::value)
  && (::pressio::is_native_container_kokkos<T2>::value
   ||  ::pressio::is_expression_acting_on_kokkos<T2>::value)
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<T, T1, T2>::value
  && (std::is_floating_point<typename ::pressio::Traits<T>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<T>::scalar_type>::value)
  && std::is_convertible<b_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<c_Type, typename ::pressio::Traits<T>::scalar_type>::value
  >
update(T & v,
    const T1 & v1, const b_Type &b,
    const T2 & v2, const c_Type &c)
{
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v1, 0));
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v2, 0));

  using scalar_t = typename ::pressio::Traits<T>::scalar_type;
  static_assert(impl::_kokkosUpdateAdmissibleOperands<scalar_t,T,T1,T2>::value,"");

  const scalar_t b_(b);
  const scalar_t c_(c);

  constexpr auto zero = static_cast<scalar_t>(0);
  if (b_ == zero) {
    ::pressio::ops::update(v, v2, c_);
  } else if (c_ == zero) {
    ::pressio::ops::update(v, v1, b_);
  } else {
    using v_t = typename impl::NativeType<T>::type;
    using v1_t = typename impl::NativeType<T1>::type;
    using v2_t = typename impl::NativeType<T2>::type;

    using fnctr_t = ::pressio::ops::impl::DoUpdateTwoTermsFunctor<v_t,v1_t,v2_t,scalar_t>;
    fnctr_t F(impl::get_native(v),
              impl::get_native(v1),
              impl::get_native(v2), b_, c_);
    Kokkos::parallel_for(v.extent(0), F);
  }
}

// //----------------------------------------------------------------------
// //  overloads for computing:
// //	V = a * V + b * V1 + c * V2 + d * V3
// //----------------------------------------------------------------------
template<
  typename T, typename T1, typename T2, typename T3,
  typename a_Type, typename b_Type, typename c_Type, typename d_Type
  >
std::enable_if_t<
  // rank-1 update common constraints
     ::pressio::Traits<T>::rank == 1
  && ::pressio::Traits<T1>::rank == 1
  && ::pressio::Traits<T2>::rank == 1
  && ::pressio::Traits<T3>::rank == 1
  // TPL/container specific
  && (::pressio::is_native_container_kokkos<T>::value
   || ::pressio::is_expression_acting_on_kokkos<T>::value)
  && (::pressio::is_native_container_kokkos<T1>::value
   ||  ::pressio::is_expression_acting_on_kokkos<T1>::value)
  && (::pressio::is_native_container_kokkos<T2>::value
   ||  ::pressio::is_expression_acting_on_kokkos<T2>::value)
  && (::pressio::is_native_container_kokkos<T3>::value
   ||  ::pressio::is_expression_acting_on_kokkos<T3>::value)
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<T, T1, T2, T3>::value
  && (std::is_floating_point<typename ::pressio::Traits<T>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<T>::scalar_type>::value)
  && std::is_convertible<a_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<b_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<c_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<d_Type, typename ::pressio::Traits<T>::scalar_type>::value
  >
update(T & v, const a_Type &a,
    const T1 & v1, const b_Type &b,
    const T2 & v2, const c_Type &c,
    const T3 & v3, const d_Type &d)
{
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v1, 0));
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v2, 0));
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v3, 0));

  using scalar_t = typename ::pressio::Traits<T>::scalar_type;
  static_assert(impl::_kokkosUpdateAdmissibleOperands<scalar_t,T,T1,T2,T3>::value,"");

  const scalar_t a_(a);
  const scalar_t b_(b);
  const scalar_t c_(c);
  const scalar_t d_(d);

  constexpr auto zero = static_cast<scalar_t>(0);
  if (b_ == zero) {
    ::pressio::ops::update(v, a_, v2, c_, v3, d_);
  } else if (c_ == zero) {
    ::pressio::ops::update(v, a_, v1, b_, v3, d_);
  } else if (d_ == zero) {
    ::pressio::ops::update(v, a_, v1, b_, v2, c_);
  } else {
    using v_t = typename impl::NativeType<T>::type;
    using v1_t = typename impl::NativeType<T1>::type;
    using v2_t = typename impl::NativeType<T2>::type;
    using v3_t = typename impl::NativeType<T3>::type;

    using fnctr_t = ::pressio::ops::impl::DoUpdateThreeTermsFunctor<v_t,v1_t,v2_t,v3_t,scalar_t>;
    fnctr_t F(impl::get_native(v),
              impl::get_native(v1),
              impl::get_native(v2),
              impl::get_native(v3),
              a_, b_, c_, d_);
    Kokkos::parallel_for(v.extent(0), F);
  }
}

template<
  typename T, typename T1, typename T2, typename T3,
  typename b_Type, typename c_Type, typename d_Type
  >
std::enable_if_t<
  // rank-1 update common constraints
     ::pressio::Traits<T>::rank == 1
  && ::pressio::Traits<T1>::rank == 1
  && ::pressio::Traits<T2>::rank == 1
  && ::pressio::Traits<T3>::rank == 1
  // TPL/container specific
  && (::pressio::is_native_container_kokkos<T>::value
   || ::pressio::is_expression_acting_on_kokkos<T>::value)
  && (::pressio::is_native_container_kokkos<T1>::value
   ||  ::pressio::is_expression_acting_on_kokkos<T1>::value)
  && (::pressio::is_native_container_kokkos<T2>::value
   ||  ::pressio::is_expression_acting_on_kokkos<T2>::value)
  && (::pressio::is_native_container_kokkos<T3>::value
   ||  ::pressio::is_expression_acting_on_kokkos<T3>::value)
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<T, T1, T2, T3>::value
  && (std::is_floating_point<typename ::pressio::Traits<T>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<T>::scalar_type>::value)
  && std::is_convertible<b_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<c_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<d_Type, typename ::pressio::Traits<T>::scalar_type>::value
  >
update(T & v,
    const T1 & v1, const b_Type &b,
    const T2 & v2, const c_Type &c,
    const T3 & v3, const d_Type &d)
{
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v1, 0));
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v2, 0));
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v3, 0));

  using scalar_t = typename ::pressio::Traits<T>::scalar_type;
  static_assert(impl::_kokkosUpdateAdmissibleOperands<scalar_t,T,T1,T2,T3>::value,"");

  const scalar_t b_(b);
  const scalar_t c_(c);
  const scalar_t d_(d);

  constexpr auto zero = static_cast<scalar_t>(0);
  if (b_ == zero) {
    ::pressio::ops::update(v, v2, c_, v3, d_);
  } else if (c_ == zero) {
    ::pressio::ops::update(v, v1, b_, v3, d_);
  } else if (d_ == zero) {
    ::pressio::ops::update(v, v1, b_, v2, c_);
  } else {
    using v_t = typename impl::NativeType<T>::type;
    using v1_t = typename impl::NativeType<T1>::type;
    using v2_t = typename impl::NativeType<T2>::type;
    using v3_t = typename impl::NativeType<T3>::type;

    using fnctr_t = ::pressio::ops::impl::DoUpdateThreeTermsFunctor<v_t,v1_t,v2_t,v3_t,scalar_t>;
    fnctr_t F(impl::get_native(v),
              impl::get_native(v1),
              impl::get_native(v2),
              impl::get_native(v3),
              b_, c_, d_);
    Kokkos::parallel_for(v.extent(0), F);
  }
}

//----------------------------------------------------------------------
//  overloads for computing:
//	V = a * V + b * V1 + c * V2 + d * V3 + e * V4
//----------------------------------------------------------------------
template<
  typename T, typename T1, typename T2, typename T3, typename T4,
  typename a_Type, typename b_Type, typename c_Type, typename d_Type, typename e_Type
  >
std::enable_if_t<
  // rank-1 update common constraints
     ::pressio::Traits<T>::rank == 1
  && ::pressio::Traits<T1>::rank == 1
  && ::pressio::Traits<T2>::rank == 1
  && ::pressio::Traits<T3>::rank == 1
  && ::pressio::Traits<T4>::rank == 1
  // TPL/container specific
  && (::pressio::is_native_container_kokkos<T>::value
   || ::pressio::is_expression_acting_on_kokkos<T>::value)
  && (::pressio::is_native_container_kokkos<T1>::value
   ||  ::pressio::is_expression_acting_on_kokkos<T1>::value)
  && (::pressio::is_native_container_kokkos<T2>::value
   ||  ::pressio::is_expression_acting_on_kokkos<T2>::value)
  && (::pressio::is_native_container_kokkos<T3>::value
   ||  ::pressio::is_expression_acting_on_kokkos<T3>::value)
  && (::pressio::is_native_container_kokkos<T4>::value
   ||  ::pressio::is_expression_acting_on_kokkos<T4>::value)
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<T, T1, T2, T3, T4>::value
  && (std::is_floating_point<typename ::pressio::Traits<T>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<T>::scalar_type>::value)
  && std::is_convertible<a_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<b_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<c_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<d_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<e_Type, typename ::pressio::Traits<T>::scalar_type>::value
  >
update(T & v, const a_Type &a,
    const T1 & v1, const b_Type &b,
    const T2 & v2, const c_Type &c,
    const T3 & v3, const d_Type &d,
    const T4 & v4, const e_Type &e)
{
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v1, 0));
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v2, 0));
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v3, 0));
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v4, 0));

  using scalar_t = typename ::pressio::Traits<T>::scalar_type;
  static_assert(impl:: _kokkosUpdateAdmissibleOperands<scalar_t,T,T1,T2,T3,T4>::value,"");

  const scalar_t a_(a);
  const scalar_t b_(b);
  const scalar_t c_(c);
  const scalar_t d_(d);
  const scalar_t e_(e);

  constexpr auto zero = static_cast<scalar_t>(0);
  if (b_ == zero) {
    ::pressio::ops::update(v, a_, v2, c_, v3, d_, v4, e_);
  } else if (c_ == zero) {
    ::pressio::ops::update(v, a_, v1, b_, v3, d_, v4, e_);
  } else if (d_ == zero) {
    ::pressio::ops::update(v, a_, v1, b_, v2, c_, v4, e_);
  } else if (e_ == zero) {
    ::pressio::ops::update(v, a_, v1, b_, v2, c_, v3, d_);
  } else {
    using v_t = typename impl::NativeType<T>::type;
    using v1_t = typename impl::NativeType<T1>::type;
    using v2_t = typename impl::NativeType<T2>::type;
    using v3_t = typename impl::NativeType<T3>::type;
    using v4_t = typename impl::NativeType<T4>::type;

    using fnctr_t = ::pressio::ops::impl::DoUpdateFourTermsFunctor<v_t,v1_t,v2_t,v3_t,v4_t,scalar_t>;
    fnctr_t F(impl::get_native(v),
              impl::get_native(v1),
              impl::get_native(v2),
              impl::get_native(v3),
              impl::get_native(v4),
              a_, b_, c_, d_, e_);
    Kokkos::parallel_for(v.extent(0), F);
  }
}

template<
  typename T, typename T1, typename T2, typename T3, typename T4,
  typename b_Type, typename c_Type, typename d_Type, typename e_Type
  >
std::enable_if_t<
  // rank-1 update common constraints
     ::pressio::Traits<T>::rank == 1
  && ::pressio::Traits<T1>::rank == 1
  && ::pressio::Traits<T2>::rank == 1
  && ::pressio::Traits<T3>::rank == 1
  && ::pressio::Traits<T4>::rank == 1
  // TPL/container specific
  && (::pressio::is_native_container_kokkos<T>::value
   || ::pressio::is_expression_acting_on_kokkos<T>::value)
  && (::pressio::is_native_container_kokkos<T1>::value
   ||  ::pressio::is_expression_acting_on_kokkos<T1>::value)
  && (::pressio::is_native_container_kokkos<T2>::value
   ||  ::pressio::is_expression_acting_on_kokkos<T2>::value)
  && (::pressio::is_native_container_kokkos<T3>::value
   ||  ::pressio::is_expression_acting_on_kokkos<T3>::value)
  && (::pressio::is_native_container_kokkos<T4>::value
   ||  ::pressio::is_expression_acting_on_kokkos<T4>::value)
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<T, T1, T2, T3, T4>::value
  && (std::is_floating_point<typename ::pressio::Traits<T>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<T>::scalar_type>::value)
  && std::is_convertible<b_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<c_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<d_Type, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<e_Type, typename ::pressio::Traits<T>::scalar_type>::value
  >
update(T & v,
    const T1 & v1, const b_Type &b,
    const T2 & v2, const c_Type &c,
    const T3 & v3, const d_Type &d,
    const T4 & v4, const e_Type &e)
{
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v1, 0));
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v2, 0));
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v3, 0));
  assert(::pressio::ops::extent(v, 0) == ::pressio::ops::extent(v4, 0));

  using scalar_t = typename ::pressio::Traits<T>::scalar_type;
  static_assert(impl::_kokkosUpdateAdmissibleOperands<scalar_t,T,T1,T2,T3,T4>::value,"");

  const scalar_t b_(b);
  const scalar_t c_(c);
  const scalar_t d_(d);
  const scalar_t e_(e);

  constexpr auto zero = static_cast<scalar_t>(0);
  if (b_ == zero) {
    ::pressio::ops::update(v, v2, c_, v3, d_, v4, e_);
  } else if (c_ == zero) {
    ::pressio::ops::update(v, v1, b_, v3, d_, v4, e_);
  } else if (d_ == zero) {
    ::pressio::ops::update(v, v1, b_, v2, c_, v4, e_);
  } else if (e_ == zero) {
    ::pressio::ops::update(v, v1, b_, v2, c_, v3, d_);
  } else {
    using v_t = typename impl::NativeType<T>::type;
    using v1_t = typename impl::NativeType<T1>::type;
    using v2_t = typename impl::NativeType<T2>::type;
    using v3_t = typename impl::NativeType<T3>::type;
    using v4_t = typename impl::NativeType<T4>::type;

    using fnctr_t = ::pressio::ops::impl::DoUpdateFourTermsFunctor<v_t,v1_t,v2_t,v3_t,v4_t,scalar_t>;
    fnctr_t F(impl::get_native(v),
              impl::get_native(v1),
              impl::get_native(v2),
              impl::get_native(v3),
              impl::get_native(v4),
              b_, c_, d_, e_);
    Kokkos::parallel_for(v.extent(0), F);
  }
}

}}//end namespace pressio::ops
#endif  // OPS_KOKKOS_OPS_VECTOR_UPDATE_HPP_
