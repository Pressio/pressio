/*
//@HEADER
// ************************************************************************
//
// solvers_legitimate_linear_solver_for_nonlinear_least_squares.hpp
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

#ifndef SOLVERS_CONSTRAINTS_SOLVERS_LEGITIMATE_LINEAR_SOLVER_FOR_NONLINEAR_LEAST_SQUARES_HPP_
#define SOLVERS_CONSTRAINTS_SOLVERS_LEGITIMATE_LINEAR_SOLVER_FOR_NONLINEAR_LEAST_SQUARES_HPP_

namespace pressio{ namespace solvers{ namespace constraints {

template <
  typename T,
  typename state_type,
  typename rhs_type = state_type,
  typename enable = void
  >
struct linear_solver_for_nonlinear_least_squares : std::false_type
{
  static_assert
  (!std::is_const<T>::value,
   "The linear solver type cannot be cv-qualified: maybe you are using a const object?");
};

template <typename T, typename state_type>
struct linear_solver_for_nonlinear_least_squares<
  T, state_type, state_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::solvers::predicates::has_matrix_typedef<T>::value and
    // the matrix_type is not void
    !std::is_void<typename T::matrix_type>::value and
    // has a solve method
    std::is_void<
      decltype
      (
       std::declval<T>().solve
       (
        std::declval<typename T::matrix_type const &>(), // A
        std::declval<state_type const &>(), // b
        std::declval<state_type &>() // x
        )
       )
      >::value
    >
  > : std::true_type{};

}}} // namespace pressio::solvers::constraints
#endif  // SOLVERS_CONSTRAINTS_SOLVERS_LEGITIMATE_LINEAR_SOLVER_FOR_NONLINEAR_LEAST_SQUARES_HPP_
