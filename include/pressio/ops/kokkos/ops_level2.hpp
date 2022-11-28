/*
//@HEADER
// ************************************************************************
//
// ops_level2.hpp
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

#ifndef OPS_KOKKOS_OPS_LEVEL2_HPP_
#define OPS_KOKKOS_OPS_LEVEL2_HPP_

#include "KokkosBlas2_gemv.hpp"

namespace pressio{ namespace ops{

/*
 * multi_vector prod vector
 *
 * y = beta * y + alpha*op(A)*x
 *
*/

//-------------------------------
// specialize for op(A) = A
//-------------------------------
template <
  class A_type, class x_type, class y_type,
  class alpha_t, class beta_t
  >
::pressio::mpl::enable_if_t<
  // level2 common constraints
  ::pressio::Traits<A_type>::rank == 2 and
  ::pressio::Traits<x_type>::rank == 1 and
  ::pressio::Traits<y_type>::rank == 1 and
  // TPL/container specific
  (::pressio::is_native_container_kokkos<A_type>::value or
   ::pressio::is_expression_acting_on_kokkos<A_type>::value) and
  (::pressio::is_native_container_kokkos<x_type>::value or
   ::pressio::is_expression_acting_on_kokkos<x_type>::value) and
  (::pressio::is_native_container_kokkos<y_type>::value or
    ::pressio::is_expression_acting_on_kokkos<y_type>::value) and
  // scalar compatibility
  ::pressio::all_have_traits_and_same_scalar<A_type, x_type, y_type>::value and
  std::is_convertible<alpha_t, typename A_type::non_const_value_type>::value and
  std::is_convertible<beta_t, typename y_type::non_const_value_type>::value and
  (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value or
   std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  >
product(::pressio::nontranspose /*unused*/,
	const alpha_t alpha,
	const A_type & A,
	const x_type & x,
	const beta_t beta,
	y_type & y)
{
  assert( y.extent(0) == A.extent(0) );
  assert( A.extent(1) == x.extent(0) );

  using sc_t = typename A_type::non_const_value_type;
  const sc_t alpha_(alpha);
  const sc_t beta_(beta);

  // 2022-11-28: Kokkos requires definition of
  // Kokkos::Details::ArithTraits<sc_t>::zero()
  // so we can use it here
  const auto zero = Kokkos::Details::ArithTraits<sc_t>::zero();
  if (static_cast<sc_t>(alpha_) == zero) {
    if (beta_ == zero) {
      ::pressio::ops::set_zero(y);
    } else {
      ::pressio::ops::scale(y, beta_);
    }
    return;
  }

  const char ctA = 'N';

  auto & y_n = impl::get_native(y);
  const auto & A_n = impl::get_native(A);
  const auto & x_n = impl::get_native(x);
  ::KokkosBlas::gemv( &ctA, alpha_, A_n, x_n, beta_, y_n);
}

//-------------------------------
// specialize for op(A) = A^T
//-------------------------------
template <
  class A_type, class x_type, class y_type,
  class alpha_t, class beta_t
  >
::pressio::mpl::enable_if_t<
  // level2 common constraints
  ::pressio::Traits<A_type>::rank == 2 and
  ::pressio::Traits<x_type>::rank == 1 and
  ::pressio::Traits<y_type>::rank == 1 and
  // TPL/container specific
  (::pressio::is_native_container_kokkos<A_type>::value or
   ::pressio::is_expression_acting_on_kokkos<A_type>::value) and
  (::pressio::is_native_container_kokkos<x_type>::value or
   ::pressio::is_expression_acting_on_kokkos<x_type>::value) and
  (::pressio::is_native_container_kokkos<y_type>::value or
    ::pressio::is_expression_acting_on_kokkos<y_type>::value) and
  // scalar compatibility
  ::pressio::all_have_traits_and_same_scalar<A_type, x_type, y_type>::value and
  std::is_convertible<alpha_t, typename A_type::non_const_value_type>::value and
  std::is_convertible<beta_t, typename y_type::non_const_value_type>::value and
  (std::is_floating_point<typename ::pressio::Traits<A_type>::scalar_type>::value or
   std::is_integral<typename ::pressio::Traits<A_type>::scalar_type>::value)
  >
product(::pressio::transpose /*unused*/,
	const alpha_t alpha,
	const A_type & A,
	const x_type & x,
	const beta_t beta,
	y_type & y)
{
  assert( y.extent(0) == A.extent(1) );
  assert( A.extent(0) == x.extent(0) );

  using sc_t = typename A_type::non_const_value_type;
  const sc_t alpha_(alpha);
  const sc_t beta_(beta);

  // 2022-11-28: Kokkos requires definition of
  // Kokkos::Details::ArithTraits<sc_t>::zero()
  // so we can use it here
  const auto zero = Kokkos::Details::ArithTraits<sc_t>::zero();
  if (static_cast<sc_t>(alpha_) == zero) {
    if (beta_ == zero) {
      ::pressio::ops::set_zero(y);
    } else {
      ::pressio::ops::scale(y, beta_);
    }
    return;
  }

  const char ctA = 'T';

  auto & y_n = impl::get_native(y);
  const auto & A_n = impl::get_native(A);
  const auto & x_n = impl::get_native(x);
  ::KokkosBlas::gemv(&ctA, alpha_, A_n, x_n, beta_, y_n);
}

}}//end namespace pressio::ops
#endif  // OPS_KOKKOS_OPS_LEVEL2_HPP_
