/*
//@HEADER
// ************************************************************************
//
// rom_has_const_apply_discrete_time_jacobian_method_accept_step_time_dt_operand_result_n_states_returning_void.hpp
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

#ifndef ROM_PREDICATES_APPLY_DISCRETE_TIME_JACOBIAN_METHODS_ROM_HAS_CONST_APPLY_DISCRETE_TIME_JACOBIAN_METHOD_ACCEPT_STEP_TIME_DT_OPERAND_RESULT_N_STATES_RETURNING_VOID_HPP_
#define ROM_PREDICATES_APPLY_DISCRETE_TIME_JACOBIAN_METHODS_ROM_HAS_CONST_APPLY_DISCRETE_TIME_JACOBIAN_METHOD_ACCEPT_STEP_TIME_DT_OPERAND_RESULT_N_STATES_RETURNING_VOID_HPP_

namespace pressio{ namespace rom{ namespace predicates {

template <
  class T, int n, class step_t, class sc_t, class state_t, class operand_t, class result_t,
  class = void
  >
struct has_const_apply_discrete_time_jacobian_method_accept_step_time_dt_operand_result_n_states_returning_void
  : std::false_type{};

template <class T, class step_t, class sc_t, class state_t, class operand_t, class result_t>
struct has_const_apply_discrete_time_jacobian_method_accept_step_time_dt_operand_result_n_states_returning_void<
  T, 1, step_t, sc_t, state_t, operand_t, result_t,
  ::pressio::mpl::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T const>().applyDiscreteTimeJacobian
       (
	std::declval<step_t const &>(),
	std::declval<sc_t const &>(),
	std::declval<sc_t const &>(),
	std::declval<operand_t const &>(),
	std::declval<result_t &>(),
	std::declval<state_t const&>()
	)
       )
      >::value
    >
  > : std::true_type{};


template <class T, class step_t, class sc_t, class state_t, class operand_t, class result_t>
struct has_const_apply_discrete_time_jacobian_method_accept_step_time_dt_operand_result_n_states_returning_void<
  T, 2, step_t, sc_t, state_t, operand_t, result_t,
  ::pressio::mpl::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T const>().applyDiscreteTimeJacobian
       (
	std::declval<step_t const &>(),
	std::declval<sc_t const &>(),
	std::declval<sc_t const &>(),
	std::declval<operand_t const &>(),
	std::declval<result_t &>(),
	std::declval<state_t const&>(),
	std::declval<state_t const&>()
	)
       )
      >::value
    >
  > : std::true_type{};


template <class T, class step_t, class sc_t, class state_t, class operand_t, class result_t>
struct has_const_apply_discrete_time_jacobian_method_accept_step_time_dt_operand_result_n_states_returning_void<
  T, 3, step_t, sc_t, state_t, operand_t, result_t,
  ::pressio::mpl::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T const>().applyDiscreteTimeJacobian
       (
	std::declval<step_t const &>(),
	std::declval<sc_t const &>(),
	std::declval<sc_t const &>(),
	std::declval<operand_t const &>(),
	std::declval<result_t &>(),
	std::declval<state_t const&>(),
	std::declval<state_t const&>(),
	std::declval<state_t const&>()
	)
       )
      >::value
    >
  > : std::true_type{};

}}}
#endif  // ROM_PREDICATES_APPLY_DISCRETE_TIME_JACOBIAN_METHODS_ROM_HAS_CONST_APPLY_DISCRETE_TIME_JACOBIAN_METHOD_ACCEPT_STEP_TIME_DT_OPERAND_RESULT_N_STATES_RETURNING_VOID_HPP_
