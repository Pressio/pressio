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

//
// rhs only
//
template<class T, class enable = void>
struct SemiDiscreteSystemWithRhs : std::false_type{};

template<class T>
struct SemiDiscreteSystemWithRhs<
  T,
  mpl::enable_if_t<
    ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_right_hand_side_typedef<T>::value
    //
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::right_hand_side_type>::value
    && ::pressio::VectorSpaceElementsWithSameField<
      typename T::state_type, typename T::right_hand_side_type
      >::value
    && std::is_convertible<
      typename T::independent_variable_type,
      typename ::pressio::Traits<typename T::state_type>::scalar_type
      >::value
    //
    && ::pressio::ode::has_const_create_state_method_return_result<
      T, typename T::state_type >::value
    && ::pressio::ode::has_const_create_rhs_method_return_result<
      T, typename T::right_hand_side_type >::value
    && ::pressio::ode::has_const_rhs_method_accept_state_indvar_result_return_void<
      T, typename T::state_type,typename T::independent_variable_type,typename T::right_hand_side_type>::value
   >
  > : std::true_type{};

//
// rhs, jacobian
//
template<class T, class enable = void>
struct SemiDiscreteSystemWithRhsAndJacobian : std::false_type{};

template<class T>
struct SemiDiscreteSystemWithRhsAndJacobian<
  T,
  mpl::enable_if_t<
    SemiDiscreteSystemWithRhs<T>::value
    && ::pressio::has_jacobian_typedef<T>::value
    && std::is_copy_constructible<typename T::jacobian_type>::value
    && ::pressio::VectorSpaceElementsWithSameField<
      typename T::state_type, typename T::jacobian_type
      >::value
    //
    // jacobian methods
    && ::pressio::ode::has_const_create_jacobian_method_return_result<
      T, typename T::jacobian_type >::value
    && ::pressio::ode::has_const_jacobian_method_accept_state_indvar_result_return_void<
      T, typename T::state_type, typename T::independent_variable_type, typename T::jacobian_type
      >::value
    >
  > : std::true_type{};

//
// rhs and varying mass matrix
//
template<class T, class enable = void>
struct SemiDiscreteSystemWithRhsAndMassMatrix : std::false_type{};

template<class T>
struct SemiDiscreteSystemWithRhsAndMassMatrix<
  T,
  mpl::enable_if_t<
    SemiDiscreteSystemWithRhs<T>::value
    && ::pressio::has_mass_matrix_typedef<T>::value
    && std::is_copy_constructible<typename T::mass_matrix_type>::value
    && ::pressio::VectorSpaceElementsWithSameField<
      typename T::state_type, typename T::mass_matrix_type
      >::value
    //
    // mass matrix
    && ::pressio::ode::has_const_create_mass_matrix_method_return_result<
      T, typename T::mass_matrix_type >::value
    && ::pressio::ode::has_const_mass_matrix_method_accept_state_indvar_result_return_void<
      T, typename T::state_type, typename T::independent_variable_type,typename T::mass_matrix_type>::value
   >
  > : std::true_type{};


//
// rhs and CONSTANT mass matrix
//
template<class T, class enable = void>
struct SemiDiscreteSystemWithRhsAndConstantMassMatrix : std::false_type{};

template<class T>
struct SemiDiscreteSystemWithRhsAndConstantMassMatrix<
  T,
  mpl::enable_if_t<
    SemiDiscreteSystemWithRhs<T>::value
    && ::pressio::has_mass_matrix_typedef<T>::value
    && std::is_copy_constructible<typename T::mass_matrix_type>::value
    && ::pressio::VectorSpaceElementsWithSameField<
      typename T::state_type, typename T::mass_matrix_type
      >::value
    //
    // mass matrix
    && ::pressio::ode::has_const_create_mass_matrix_method_return_result<
      T, typename T::mass_matrix_type >::value
    && ::pressio::ode::has_const_mass_matrix_method_accept_result_return_void<
      T, typename T::mass_matrix_type>::value
   >
  > : std::true_type{};


//
// rhs, jac, mass matrix
//
template<class T, class enable = void>
struct CompleteSemiDiscreteSystem : std::false_type{};

template<class T>
struct CompleteSemiDiscreteSystem<
  T,
  mpl::enable_if_t<
    SemiDiscreteSystemWithRhsAndJacobian<T>::value
    && SemiDiscreteSystemWithRhsAndMassMatrix<T>::value
    >
  > : std::true_type{};

//
// rhs, jac, CONSTANT mass matrix
//
template<class T, class enable = void>
struct CompleteSemiDiscreteSystemWithConstantMassMatrix : std::false_type{};

template<class T>
struct CompleteSemiDiscreteSystemWithConstantMassMatrix<
  T,
  mpl::enable_if_t<
    SemiDiscreteSystemWithRhsAndJacobian<T>::value
    && SemiDiscreteSystemWithRhsAndConstantMassMatrix<T>::value
    >
  > : std::true_type{};

//
// fully discrete
//
template<class T, int NumStates, class enable = void>
struct FullyDiscreteSystemWithJacobian : std::false_type{};

template<class T, int NumStates>
struct FullyDiscreteSystemWithJacobian<
  T, NumStates,
  mpl::enable_if_t<
       ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_discrete_residual_typedef<T>::value
    && ::pressio::has_discrete_jacobian_typedef<T>::value
    //
    && std::is_copy_constructible<typename T::discrete_residual_type>::value
    && std::is_copy_constructible<typename T::discrete_jacobian_type>::value
    && ::pressio::VectorSpaceElementsWithSameField<
	 typename T::state_type, typename T::discrete_residual_type, typename T::discrete_jacobian_type
	 >::value
    && std::is_convertible<
	 typename T::independent_variable_type,
	 typename ::pressio::Traits<typename T::state_type>::scalar_type
	 >::value
    //
    // creation methods
    && ::pressio::ode::has_const_create_state_method_return_result<
      T, typename T::state_type >::value
    && ::pressio::ode::has_const_create_discrete_residual_method_return_result<
      T, typename T::discrete_residual_type>::value
    && ::pressio::ode::has_const_create_discrete_jacobian_method_return_result<
      T, typename T::discrete_jacobian_type>::value
    //
    // discrete residual
    && ::pressio::ode::has_const_discrete_residual_method_accept_step_indvar_dt_result_n_states_return_void<
      T, NumStates, int, typename T::independent_variable_type, typename T::state_type,
      typename T::discrete_residual_type>::value
    //
    // discrete jacobian
    && ::pressio::ode::has_const_discrete_jacobian_method_accept_step_indvar_dt_result_n_states_return_void<
      T, NumStates, int, typename T::independent_variable_type, typename T::state_type,
      typename T::discrete_jacobian_type>::value
    >
  > : std::true_type{};

}} // namespace pressio::ode::constraints
#endif  // ODE_STEPPERS_CONSTRAINTS_ODE_SYSTEM_API_HPP_
