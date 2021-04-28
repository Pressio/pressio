/*
//@HEADER
// ************************************************************************
//
// solvers_ops_normal_equations_rj_api.hpp
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

#ifndef SOLVERS_CONSTRAINTS_SOLVERS_OPS_NORMAL_EQUATIONS_RJ_API_HPP_
#define SOLVERS_CONSTRAINTS_SOLVERS_OPS_NORMAL_EQUATIONS_RJ_API_HPP_

namespace pressio{ namespace solvers{ namespace constraints {

template <
  typename T,
  typename sc_t,
  typename H_t,
  typename g_t,
  typename J_t,
  typename r_t,
  typename enable = void
>
struct ops_normal_equations_rj_api : std::false_type{};

template <
  typename T,
  typename sc_t,
  typename H_t,
  typename g_t,
  typename J_t,
  typename r_t
>
struct ops_normal_equations_rj_api<
  T, sc_t, H_t, g_t, J_t, r_t,
  ::pressio::mpl::enable_if_t<

    // 1. need non-void product returning J^T J
    std::is_same<
      decltype
      (
       std::declval< T const &>().template product<H_t>
       (
        std::declval< ::pressio::transpose >(),
        std::declval< ::pressio::nontranspose >(),
        std::declval< sc_t>(),
        std::declval<
          typename containers::details::traits<J_t>::wrapped_t const & >(),
        std::declval<
          typename containers::details::traits<J_t>::wrapped_t const & >()
        )
       ),
      H_t
      >::value
    and

    // 2. need void product computing  J^T J
    std::is_void<
      decltype
      (
       std::declval< T const &>().product
       (
        std::declval< ::pressio::transpose >(),
        std::declval< ::pressio::nontranspose >(),
        std::declval< sc_t>(),
        std::declval< typename containers::details::traits<J_t>::wrapped_t const & >(),
        std::declval< typename containers::details::traits<J_t>::wrapped_t const & >(),
        std::declval< sc_t>(),
        std::declval< H_t & >()
        )
       )
      >::value
    and

    // 3. need void product computing J^T r
    std::is_void<
      decltype
      (
       std::declval< T const &>().product
       (
        std::declval< ::pressio::transpose >(),
        std::declval< sc_t>(),
        std::declval< typename containers::details::traits<J_t>::wrapped_t const & >(),
        std::declval< typename containers::details::traits<r_t>::wrapped_t const & >(),
        std::declval< sc_t>(),
        std::declval< g_t & >()
        )
       )
      >::value
    and

    // 4. need norm2 method
    std::is_same<
      decltype
      (
       std::declval< T const &>().norm2
       (
        std::declval< typename containers::details::traits<r_t>::wrapped_t const & >()
        )
       ), sc_t
      >::value
   >
  > : std::true_type{};

}}} // namespace pressio::solvers::constraints
#endif  // SOLVERS_CONSTRAINTS_SOLVERS_OPS_NORMAL_EQUATIONS_RJ_API_HPP_
