/*
//@HEADER
// ************************************************************************
//
// ops_elementwise_multiply.hpp
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

#ifndef OPS_KOKKOS_OPS_ELEMENTWISE_MULTIPLY_HPP_
#define OPS_KOKKOS_OPS_ELEMENTWISE_MULTIPLY_HPP_

#include <KokkosBlas1_mult.hpp>

namespace pressio{ namespace ops{

//----------------------------------------------------------------------
// computing elementwise:  y = beta * y + alpha * x * z
//----------------------------------------------------------------------
template <typename T, typename T1, typename T2, typename alpha_t, typename beta_t>
::pressio::mpl::enable_if_t<
  // common elementwise_multiply constraints
     ::pressio::Traits<T>::rank == 1
  && ::pressio::Traits<T1>::rank == 1
  && ::pressio::Traits<T2>::rank == 1
  // TPL/container specific
  && (::pressio::is_native_container_kokkos<T>::value
   || ::pressio::is_expression_acting_on_kokkos<T>::value)
  && (::pressio::is_native_container_kokkos<T1>::value
   || ::pressio::is_expression_acting_on_kokkos<T1>::value)
  && (::pressio::is_native_container_kokkos<T2>::value
   || ::pressio::is_expression_acting_on_kokkos<T2>::value)
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<T, T1, T2>::value
  && std::is_convertible<alpha_t, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<beta_t,  typename ::pressio::Traits<T>::scalar_type>::value
  && (std::is_floating_point<typename ::pressio::Traits<T>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<T>::scalar_type>::value)
  >
elementwise_multiply(alpha_t alpha,
		     const T & x,
		     const T1 & z,
		     beta_t beta,
		     T2 & y)
{
  assert(x.extent(0) == z.extent(0));
  assert(z.extent(0) == y.extent(0));

  using sc_t = typename ::pressio::Traits<T>::scalar_type;
  const sc_t alpha_(alpha);
  const sc_t beta_(beta);

  KokkosBlas::mult(beta_, impl::get_native(y),
    alpha_, impl::get_native(x), impl::get_native(z) );
}

}}//end namespace pressio::ops
#endif  // OPS_KOKKOS_OPS_ELEMENTWISE_MULTIPLY_HPP_
