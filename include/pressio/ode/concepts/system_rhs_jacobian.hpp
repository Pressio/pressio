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

#ifndef ODE_STEPPERS_SYSTEM_RHS_AND_JACOBIAN_HPP_
#define ODE_STEPPERS_SYSTEM_RHS_AND_JACOBIAN_HPP_

namespace pressio{ namespace ode{

#ifdef PRESSIO_ENABLE_CXX20
template <class T>
concept SystemWithRhsAndJacobian =
  std::copy_constructible<T>
  /*
    required nested aliases
  */
  && requires(){
    typename T::independent_variable_type;
    typename T::state_type;
    typename T::right_hand_side_type;
    typename T::jacobian_type;
  }
  /*
    requirements on the nested aliases
  */
  && ::pressio::ops::is_known_data_type<typename T::state_type>::value
  && ::pressio::ops::is_known_data_type<typename T::right_hand_side_type>::value
  && ::pressio::ops::is_known_data_type<typename T::jacobian_type>::value
  && all_have_traits_and_same_scalar<
	typename T::state_type,
	typename T::right_hand_side_type,
	typename T::jacobian_type
     >::value
  && std::convertible_to<
	typename T::independent_variable_type,
	scalar_trait_t<typename T::state_type>>
  /*
    compound requirements
  */
  && requires(const T & A,
	      const typename T::state_type & state,
	      const typename T::independent_variable_type & evalValue,
	      typename T::right_hand_side_type & rhs,
	      typename T::jacobian_type & jac,
	      bool computeJacobian)
  {
    { A.createState()         } -> std::same_as<typename T::state_type>;
    { A.createRightHandSide() } -> std::same_as<typename T::right_hand_side_type>;
    { A.createJacobian()      } -> std::same_as<typename T::jacobian_type>;
    { A(state, evalValue, rhs,
	jac, computeJacobian) } -> std::same_as<void>;
  };
#endif //PRESSIO_ENABLE_CXX20

}} // end namespace pressio::ode


#if not defined PRESSIO_ENABLE_CXX20

#include "./predicates/ode_has_const_create_state_method_return_result.hpp"
#include "./predicates/ode_has_const_create_rhs_method_return_result.hpp"
#include "./predicates/ode_has_const_create_jacobian_method_return_result.hpp"

namespace pressio{ namespace ode{

template<class T, class enable = void>
struct SystemWithRhsAndJacobian : std::false_type{};

template<class T>
struct SystemWithRhsAndJacobian<
  T,
  mpl::enable_if_t<
    std::is_copy_constructible<T>::value
    //
    && ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_right_hand_side_typedef<T>::value
    && ::pressio::has_jacobian_typedef<T>::value
    //
    && ::pressio::ops::is_known_data_type<typename T::state_type>::value
    && ::pressio::ops::is_known_data_type<typename T::right_hand_side_type>::value
    && ::pressio::ops::is_known_data_type<typename T::jacobian_type>::value
    && all_have_traits_and_same_scalar<
	 typename T::state_type,
	 typename T::right_hand_side_type,
	 typename T::jacobian_type>::value
    && std::is_convertible<
	 typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type>>::value
    //
    && ::pressio::ode::has_const_create_state_method_return_result<
      T, typename T::state_type >::value
    && ::pressio::ode::has_const_create_rhs_method_return_result<
      T, typename T::right_hand_side_type >::value
    && ::pressio::ode::has_const_create_jacobian_method_return_result<
      T, typename T::jacobian_type >::value
    //
    && std::is_void<
      decltype(
	       std::declval<T const>()
	       (
		std::declval<typename T::state_type const&>(),
		std::declval<typename T::independent_variable_type const &>(),
		std::declval<typename T::right_hand_side_type &>(),
		std::declval<typename T::jacobian_type &>(),
		std::declval<bool>()
	       )
	   )
      >::value
   >
  > : std::true_type{};

}} // end namespace pressio::ode
#endif

#endif
