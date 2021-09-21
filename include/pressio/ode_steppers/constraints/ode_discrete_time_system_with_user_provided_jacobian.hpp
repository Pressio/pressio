/*
//@HEADER
// ************************************************************************
//
// ode_discrete_time_system_with_user_provided_jacobian.hpp
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


#ifndef ODE_STEPPERS_CONSTRAINTS_ODE_DISCRETE_TIME_SYSTEM_WITH_USER_PROVIDED_JACOBIAN_HPP_
#define ODE_STEPPERS_CONSTRAINTS_ODE_DISCRETE_TIME_SYSTEM_WITH_USER_PROVIDED_JACOBIAN_HPP_

namespace pressio{ namespace ode{

template<class T, int num_states, class enable = void>
struct discrete_time_system_with_user_provided_jacobian : std::false_type{};

template<class T, int num_states>
struct discrete_time_system_with_user_provided_jacobian<
  T, num_states,
  mpl::enable_if_t<
    ::pressio::has_scalar_typedef<T>::value and
    ::pressio::has_state_typedef<T>::value and
    ::pressio::has_discrete_time_residual_typedef<T>::value and
    ::pressio::has_discrete_time_jacobian_typedef<T>::value and
    //
    // time-discrete residual
    ::pressio::ode::has_const_create_discrete_time_residual_method_return_result<
        T, typename T::discrete_time_residual_type>::value and
    ::pressio::ode::has_const_discrete_time_residual_method_accept_step_time_dt_result_n_states_return_void<
        T, num_states, int, typename T::scalar_type, typename T::state_type,
        typename T::discrete_time_residual_type>::value and
    // time-discrete jacobian
    ::pressio::ode::has_const_create_discrete_time_jacobian_method_return_result<
        T, typename T::discrete_time_jacobian_type>::value and
    ::pressio::ode::has_const_discrete_time_jacobian_method_accept_step_time_dt_result_n_states_return_void<
        T, num_states, int, typename T::scalar_type, typename T::state_type,
        typename T::discrete_time_jacobian_type>::value
    >
  > : std::true_type{};

}} // namespace pressio::ode::constraints
#endif  // ODE_STEPPERS_CONSTRAINTS_ODE_DISCRETE_TIME_SYSTEM_WITH_USER_PROVIDED_JACOBIAN_HPP_
