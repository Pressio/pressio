/*
//@HEADER
// ************************************************************************
//
// rom_discrete_time_system_with_user_provided_apply_jacobian.hpp
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

#ifndef ROM_CONSTRAINTS_SYSTEM_ROM_DISCRETE_TIME_SYSTEM_WITH_USER_PROVIDED_APPLY_JACOBIAN_HPP_
#define ROM_CONSTRAINTS_SYSTEM_ROM_DISCRETE_TIME_SYSTEM_WITH_USER_PROVIDED_APPLY_JACOBIAN_HPP_

namespace pressio{ namespace rom{ namespace constraints {

template<typename T, typename apply_jac_operand_t, typename enable = void>
struct discrete_time_system_with_user_provided_apply_jacobian : std::false_type{};

template<typename T, typename apply_jac_operand_t>
struct discrete_time_system_with_user_provided_apply_jacobian<
  T, apply_jac_operand_t,
  mpl::enable_if_t<
    !::pressio::containers::predicates::is_wrapper<apply_jac_operand_t>::value and
    ::pressio::containers::predicates::has_scalar_typedef<T>::value and
    ::pressio::ode::predicates::has_state_typedef<T>::value and
    ::pressio::ode::predicates::has_discrete_time_residual_typedef<T>::value and
    //
    ::pressio::ode::predicates::has_const_create_discrete_time_residual_method_return_result<
      T, typename T::discrete_time_residual_type
      >::value and
    ::pressio::ode::predicates::has_const_discrete_time_residual_method_accept_step_time_dt_result_n_states_return_void<
      T, 2,
      ::pressio::ode::types::step_t,
      typename T::scalar_type, typename T::state_type, typename T::discrete_time_residual_type
      >::value and
    ::pressio::ode::predicates::has_const_discrete_time_residual_method_accept_step_time_dt_result_n_states_return_void<
      T, 3, ::pressio::ode::types::step_t, typename T::scalar_type,
      typename T::state_type, typename T::discrete_time_residual_type
      >::value and
    //
    ::pressio::rom::predicates::has_const_create_apply_discrete_time_jacobian_result_method_accept_operand_return_result<
      T, apply_jac_operand_t
      >::value and
    ::pressio::rom::predicates::has_const_apply_discrete_time_jacobian_method_accept_step_time_dt_operand_result_n_states_returning_void<
     T, 2, ::pressio::ode::types::step_t, typename T::scalar_type,
     typename T::state_type, apply_jac_operand_t, apply_jac_operand_t
     >::value and
  ::pressio::rom::predicates::has_const_apply_discrete_time_jacobian_method_accept_step_time_dt_operand_result_n_states_returning_void<
    T, 3, ::pressio::ode::types::step_t,
    typename T::scalar_type, typename T::state_type, apply_jac_operand_t, apply_jac_operand_t
    >::value
  >
  > : std::true_type{};


template<typename T, typename apply_jac_operand_t>
struct discrete_time_system_with_user_provided_apply_jacobian<
  T, apply_jac_operand_t,
  mpl::enable_if_t<
    ::pressio::containers::predicates::is_wrapper<apply_jac_operand_t>::value and
    discrete_time_system_with_user_provided_apply_jacobian<
      T, typename apply_jac_operand_t::traits::wrapped_t
      >::value
    >
  > : std::true_type{};

} // namespace pressio::rom::constraints

template <typename T, typename apply_jac_operand_t>
struct why_not_discrete_time_system_with_user_provided_apply_jacobian
{
  static_assert
    (::pressio::containers::predicates::is_wrapper<apply_jac_operand_t>::value,
     "the why not checks need apply_jac_operand_t to be a wrapper type");

  static_assert
    (::pressio::containers::predicates::has_scalar_typedef<T>::value,
     "Your discrete-time adapter class is without (or has a wrong) scalar typedef");

  static_assert
    (::pressio::ode::predicates::has_state_typedef<T>::value,
     "Your discrete-time adapter class is without (or has a wrong) state typedef");

  static_assert
    (::pressio::ode::predicates::has_discrete_time_residual_typedef<T>::value,
     "Your discrete-time adapter class is without (or has a wrong) discrete time residual typedef");

  static_assert
    (::pressio::ode::predicates::has_const_create_discrete_time_residual_method_return_result<
     T, typename T::discrete_time_residual_type>::value,
     "Your discrete-time adapter class is without (or has a wrong) create discrete time residual method");

  static_assert
    (::pressio::rom::predicates::has_const_create_apply_discrete_time_jacobian_result_method_accept_operand_return_result<
     T,  typename apply_jac_operand_t::traits::wrapped_t>::value,
     "Your discrete-time adapter class is without (or has a wrong) create apply discrete time jacobian result method");

  static_assert
    (::pressio::ode::predicates::has_const_discrete_time_residual_method_accept_step_time_dt_result_n_states_return_void<
     T, 2, ::pressio::ode::types::step_t,
     typename T::scalar_type, typename T::state_type, typename T::discrete_time_residual_type>::value,
     "Your discrete-time adapter class is without (or has a wrong) discrete time residual method accepting 2 state");

  static_assert
    (::pressio::ode::predicates::has_const_discrete_time_residual_method_accept_step_time_dt_result_n_states_return_void<
     T, 3, ::pressio::ode::types::step_t,
     typename T::scalar_type, typename T::state_type, typename T::discrete_time_residual_type>::value,
     "Your discrete-time adapter class is without (or has a wrong) discrete time residual method accepting 3 states");

  static_assert
    (::pressio::rom::predicates::has_const_apply_discrete_time_jacobian_method_accept_step_time_dt_operand_result_n_states_returning_void<
     T, 2,
     ::pressio::ode::types::step_t,
     typename T::scalar_type,
     typename T::state_type,
     typename apply_jac_operand_t::traits::wrapped_t,
     typename apply_jac_operand_t::traits::wrapped_t
     >::value,
     "Your discrete-time adapter class is without (or has a wrong) apply discrete time jacobian method accepting 2 states");

  static_assert
    (::pressio::rom::predicates::has_const_apply_discrete_time_jacobian_method_accept_step_time_dt_operand_result_n_states_returning_void<
     T, 3,
     ::pressio::ode::types::step_t,
     typename T::scalar_type,
     typename T::state_type,
     typename apply_jac_operand_t::traits::wrapped_t,
     typename apply_jac_operand_t::traits::wrapped_t
     >::value,
     "Your discrete-time adapter class is without (or has a wrong) apply discrete time jacobian result method accepting 3 states");

  static constexpr bool value = true;
};

}} // namespace pressio::rom
#endif  // ROM_CONSTRAINTS_SYSTEM_ROM_DISCRETE_TIME_SYSTEM_WITH_USER_PROVIDED_APPLY_JACOBIAN_HPP_
