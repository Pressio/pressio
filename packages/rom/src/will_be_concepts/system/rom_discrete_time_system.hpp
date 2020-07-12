/*
//@HEADER
// ************************************************************************
//
// rom_discrete_time_system.hpp
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

#ifndef rom_discrete_time_system_HPP_
#define rom_discrete_time_system_HPP_

namespace pressio{ namespace rom{ namespace concepts {

template<typename T, typename enable = void>
struct discrete_time_system : std::false_type{};

// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
// template<typename T>
// struct rom_discrete_time_system<
//   T,
//   mpl::enable_if_t<
//     ::pressio::mpl::is_same<T, pybind11::object>::value
//     >
//   > : std::true_type{};
// #endif

template<typename T>
struct discrete_time_system<
  T,
  mpl::enable_if_t<
    ::pressio::containers::predicates::has_scalar_typedef<T>::value and
    ::pressio::ode::predicates::has_state_typedef<T>::value and
    ::pressio::ode::predicates::has_discrete_time_residual_typedef<T>::value and
    ::pressio::rom::predicates::has_dense_matrix_typedef<T>::value 
    and
    ///////////////////////////
    // time-discrete residual 
    ::pressio::ode::predicates::has_const_create_discrete_time_residual_method_return_result<
        T, typename T::discrete_time_residual_type>::value 
    and 
    ::pressio::ode::predicates::has_const_discrete_time_residual_method_accept_step_time_dt_result_norm_n_states_return_void<
        T, 1, ::pressio::ode::types::step_t, 
        typename T::scalar_type, 
        typename T::state_type, 
        typename T::discrete_time_residual_type>::value 
    and 
    ::pressio::ode::predicates::has_const_discrete_time_residual_method_accept_step_time_dt_result_norm_n_states_return_void<
        T, 2, ::pressio::ode::types::step_t, 
        typename T::scalar_type, 
        typename T::state_type, 
        typename T::discrete_time_residual_type>::value 
    and
    // 
    ///////////////////////////
    // apply time-discrete jacobian
    ::pressio::rom::predicates::has_const_create_apply_discrete_time_jacobian_result_method_accept_operand_return_result<
        T,  typename T::dense_matrix_type, typename T::dense_matrix_type>::value 
    and 
    ::pressio::rom::predicates::has_const_apply_discrete_time_jacobian_method_accept_step_time_dt_operand_result_n_states_returning_void<
        T, 1, 
        ::pressio::ode::types::step_t, 
        typename T::scalar_type, 
        typename T::state_type, 
        typename T::dense_matrix_type,
        typename T::dense_matrix_type>::value 
    and 
    ::pressio::rom::predicates::has_const_apply_discrete_time_jacobian_method_accept_step_time_dt_operand_result_n_states_returning_void<
        T, 2, 
        ::pressio::ode::types::step_t, 
        typename T::scalar_type, 
        typename T::state_type, 
        typename T::dense_matrix_type,
        typename T::dense_matrix_type>::value 
    >
  > : std::true_type{};


}}} // namespace pressio::rom::concepts
#endif