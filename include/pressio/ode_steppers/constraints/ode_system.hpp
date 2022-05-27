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


#ifndef ODE_STEPPERS_CONSTRAINTS_ODE_SYSTEM_API_HPP_
#define ODE_STEPPERS_CONSTRAINTS_ODE_SYSTEM_API_HPP_

namespace pressio{ namespace ode{

template<class T, class enable = void>
struct system_has_complete_mass_matrix_api : std::false_type{};

template<class T>
struct system_has_complete_mass_matrix_api<
  T,
  mpl::enable_if_t<
    ::pressio::has_time_typedef<T>::value
    and ::pressio::has_state_typedef<T>::value
    and ::pressio::has_mass_matrix_typedef<T>::value
    and ::pressio::ode::has_const_create_mass_matrix_method_return_result<
      T, typename T::mass_matrix_type >::value
    and ::pressio::ode::has_const_mass_matrix_method_accept_state_time_result_return_void<
      T, typename T::state_type, typename T::time_type,typename T::mass_matrix_type>::value
   >
  > : std::true_type{};


//
// (1) continuous, velocity only
//
template<class T, class enable = void>
struct continuous_time_system_with_at_least_velocity : std::false_type{};

template<class T>
struct continuous_time_system_with_at_least_velocity<
  T,
  mpl::enable_if_t<
    ::pressio::has_time_typedef<T>::value
    and ::pressio::has_state_typedef<T>::value
    and ::pressio::has_velocity_typedef<T>::value
    //
    and ::pressio::ode::has_const_create_state_method_return_result<
      T, typename T::state_type >::value
    //
    // velocity
    and ::pressio::ode::has_const_create_velocity_method_return_result<
      T, typename T::velocity_type >::value
    and ::pressio::ode::has_const_velocity_method_accept_state_time_result_return_void<
      T, typename T::state_type,typename T::time_type,typename T::velocity_type>::value
   >
  > : std::true_type{};

//
// (2) continuous, velocity, mass matrix
//
template<class T, class enable = void>
struct continuous_time_system_with_mass_matrix_and_at_least_velocity : std::false_type{};

template<class T>
struct continuous_time_system_with_mass_matrix_and_at_least_velocity<
  T,
  mpl::enable_if_t<
    continuous_time_system_with_at_least_velocity<T>::value
    and ::pressio::has_mass_matrix_typedef<T>::value
    //
    // mass matrix
    and ::pressio::ode::has_const_create_mass_matrix_method_return_result<
      T, typename T::mass_matrix_type >::value
    and ::pressio::ode::has_const_mass_matrix_method_accept_state_time_result_return_void<
      T, typename T::state_type, typename T::time_type,typename T::mass_matrix_type>::value
   >
  > : std::true_type{};


//
// (3) continuous, velocity, jacobian
//
template<class T, class enable = void>
struct continuous_time_system_with_user_provided_jacobian : std::false_type{};

template<class T>
struct continuous_time_system_with_user_provided_jacobian<
  T,
  mpl::enable_if_t<
    ::pressio::has_time_typedef<T>::value
    and ::pressio::has_state_typedef<T>::value
    and ::pressio::has_velocity_typedef<T>::value
    and ::pressio::has_jacobian_typedef<T>::value
    //
    and ::pressio::ode::has_const_create_state_method_return_result<
      T, typename T::state_type >::value
    //
    // residual methods
    and ::pressio::ode::has_const_create_velocity_method_return_result<
      T, typename T::velocity_type >::value
    and ::pressio::ode::has_const_velocity_method_accept_state_time_result_return_void<
      T, typename T::state_type, typename T::time_type, typename T::velocity_type
      >::value
    //
    // jacobian methods
    and ::pressio::ode::has_const_create_jacobian_method_return_result<
      T, typename T::jacobian_type >::value
    and ::pressio::ode::has_const_jacobian_method_accept_state_time_result_return_void<
      T, typename T::state_type, typename T::time_type, typename T::jacobian_type
      >::value
    >
  > : std::true_type{};

//
// (4) continuous, velocity, jac, mass matrix
//
template<class T, class enable = void>
struct continuous_time_system_with_mass_matrix_and_user_provided_jacobian : std::false_type{};

template<class T>
struct continuous_time_system_with_mass_matrix_and_user_provided_jacobian<
  T,
  mpl::enable_if_t<
    continuous_time_system_with_user_provided_jacobian<T>::value
    and ::pressio::has_mass_matrix_typedef<T>::value
    //
    // mass matrix
    and ::pressio::ode::has_const_create_mass_matrix_method_return_result<
      T, typename T::mass_matrix_type >::value
    and ::pressio::ode::has_const_mass_matrix_method_accept_state_time_result_return_void<
      T, typename T::state_type,typename T::time_type,typename T::mass_matrix_type>::value
    >
  > : std::true_type{};


//
// discrete
//
template<class T, int NumStates, class enable = void>
struct discrete_time_system_with_user_provided_jacobian : std::false_type{};

template<class T, int NumStates>
struct discrete_time_system_with_user_provided_jacobian<
  T, NumStates,
  mpl::enable_if_t<
    ::pressio::has_time_typedef<T>::value
    and ::pressio::has_state_typedef<T>::value
    and ::pressio::has_discrete_time_residual_typedef<T>::value
    and ::pressio::has_discrete_time_jacobian_typedef<T>::value
    //
    and ::pressio::ode::has_const_create_state_method_return_result<
      T, typename T::state_type >::value
    //
    // time-discrete residual
    and ::pressio::ode::has_const_create_discrete_time_residual_method_return_result<
      T, typename T::discrete_time_residual_type>::value
    and ::pressio::ode::has_const_discrete_time_residual_method_accept_step_time_dt_result_n_states_return_void<
      T, NumStates, int, typename T::time_type, typename T::state_type,
      typename T::discrete_time_residual_type>::value
    //
    // time-discrete jacobian
    and ::pressio::ode::has_const_create_discrete_time_jacobian_method_return_result<
      T, typename T::discrete_time_jacobian_type>::value
    and ::pressio::ode::has_const_discrete_time_jacobian_method_accept_step_time_dt_result_n_states_return_void<
      T, NumStates, int, typename T::time_type, typename T::state_type,
      typename T::discrete_time_jacobian_type>::value
    >
  > : std::true_type{};

}} // namespace pressio::ode::constraints
#endif  // ODE_STEPPERS_CONSTRAINTS_ODE_SYSTEM_API_HPP_
