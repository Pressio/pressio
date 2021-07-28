/*
//@HEADER
// ************************************************************************
//
// solvers_least_squares_weighting_operator.hpp
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

#ifndef SOLVERS_CONSTRAINTS_SOLVERS_LEAST_SQUARES_WEIGHTING_OPERATOR_HPP_
#define SOLVERS_CONSTRAINTS_SOLVERS_LEAST_SQUARES_WEIGHTING_OPERATOR_HPP_

namespace pressio{ namespace solvers{ namespace constraints {

template <
  typename T,
  typename r_t,
  typename j_t,
  typename enable = void
>
struct least_squares_weighting_operator_accepting_wrappers : std::false_type{};

template <typename T, typename r_t, typename j_t>
struct least_squares_weighting_operator_accepting_wrappers<
  T, r_t, j_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::predicates::are_wrappers<r_t,j_t>::value
    and
    std::is_void<
      decltype
      (
       std::declval<T const>()
       (
	std::declval<r_t const &>(),
	std::declval<r_t &>()
	)
       )
      >::value
    and
    std::is_void<
      decltype
      (
       std::declval<T const>()
       (
	std::declval<j_t const &>(),
	std::declval<j_t &>()
	)
       )
      >::value
    >
  > : std::true_type{};



template <
  typename T,
  typename r_t,
  typename j_t,
  typename enable = void
>
struct least_squares_weighting_operator_accepting_native : std::false_type{};

template <typename T, typename r_t, typename j_t>
struct least_squares_weighting_operator_accepting_native<
  T, r_t, j_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::predicates::are_wrappers<r_t,j_t>::value
    and
    std::is_void<
      decltype
      (
       std::declval<T const>()
       (
	std::declval<typename ::pressio::containers::details::traits<r_t>::wrapped_t const &>(),
	std::declval<typename ::pressio::containers::details::traits<r_t>::wrapped_t &>()
	)
       )
      >::value
    and
    std::is_void<
      decltype
      (
       std::declval<T const>()
       (
	std::declval<typename ::pressio::containers::details::traits<j_t>::wrapped_t const &>(),
	std::declval<typename ::pressio::containers::details::traits<j_t>::wrapped_t &>()
	)
       )
      >::value
    >
  > : std::true_type{};



template <
  typename T,
  typename r_t,
  typename j_t,
  typename enable = void
>
struct least_squares_weighting_operator : std::false_type{};

// specialize for when the operator takes wrappers
template <typename T, typename r_t, typename j_t>
struct least_squares_weighting_operator<
  T, r_t, j_t,
  ::pressio::mpl::enable_if_t<
  least_squares_weighting_operator_accepting_wrappers<T, r_t, j_t>::value or
  least_squares_weighting_operator_accepting_native<T, r_t, j_t>::value
  >
  > : std::true_type{};


}}} // namespace pressio::solvers::constraints
#endif  // SOLVERS_CONSTRAINTS_SOLVERS_LEAST_SQUARES_WEIGHTING_OPERATOR_HPP_
