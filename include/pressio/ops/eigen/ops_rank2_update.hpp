/*
//@HEADER
// ************************************************************************
//
// ops_rank2_update.hpp
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

#ifndef OPS_EIGEN_OPS_RANK2_UPDATE_HPP_
#define OPS_EIGEN_OPS_RANK2_UPDATE_HPP_

namespace pressio{ namespace ops{

/*
   below we constrain in all cases via is_convertible
   because the implementations are using Eigen native operations
   which are based on expressions and require
   coefficients to be convertible to scalar types of the vectors operands
 */

//----------------------------------------------------------------------
// M = a * M + b * M1
//----------------------------------------------------------------------
template<typename T, typename T1, class alpha_t, class beta_t>
std::enable_if_t<
  // rank-1 update common constraints
     ::pressio::Traits<T>::rank == 2
  && ::pressio::Traits<T1>::rank == 2
  // TPL/container specific
  && (::pressio::is_native_container_eigen<T>::value
   || ::pressio::is_expression_acting_on_eigen<T>::value)
  && (::pressio::is_native_container_eigen<T1>::value
   || ::pressio::is_expression_acting_on_eigen<T1>::value)
  // scalar compatibility
  && ::pressio::all_have_traits_and_same_scalar<T, T1>::value
  && (std::is_floating_point<typename ::pressio::Traits<T>::scalar_type>::value
   || std::is_integral<typename ::pressio::Traits<T>::scalar_type>::value)
  && std::is_convertible<alpha_t, typename ::pressio::Traits<T>::scalar_type>::value
  && std::is_convertible<beta_t, typename ::pressio::Traits<T>::scalar_type>::value
  >
update(T & M,         const alpha_t & a,
       const T1 & M1, const beta_t & b)
{
  assert(::pressio::ops::extent(M, 0) == ::pressio::ops::extent(M1, 0));
  assert(::pressio::ops::extent(M, 1) == ::pressio::ops::extent(M1, 1));

  using sc_t = typename ::pressio::Traits<T>::scalar_type;
  const sc_t a_(a);
  const sc_t b_(b);

  const auto zero = ::pressio::utils::Constants<sc_t>::zero();
  if (b_ == zero) {
    ::pressio::ops::scale(M, a_);
    return;
  }

  auto & M_n = impl::get_native(M);
  const auto & M_n1 = impl::get_native(M1);
  if (a_ == zero) {
    M_n = b_*M_n1;
  } else {
    M_n = a_*M_n + b_*M_n1;
  }
}

}}//end namespace pressio::ops
#endif  // OPS_EIGEN_OPS_RANK2_UPDATE_HPP_
