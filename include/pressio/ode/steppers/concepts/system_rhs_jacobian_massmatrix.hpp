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

#ifndef ODE_STEPPERS_SYSTEM_RHS_AND_JACOBIAN_AND_MASS_MATRIX_HPP_
#define ODE_STEPPERS_SYSTEM_RHS_AND_JACOBIAN_AND_MASS_MATRIX_HPP_

#ifdef PRESSIO_ENABLE_CXX20
#include <concepts>


// this is here so that we can clearly show it in the
// doc via rst literal include directive
namespace pressio{ namespace ode{

template <class T>
concept SystemWithRhsJacobianMassMatrix =
  std::copy_constructible<T>
  /* must have nested aliases */
  && requires(){
    typename T::independent_variable_type;
    typename T::state_type;
    typename T::right_hand_side_type;
    typename T::jacobian_type;
    typename T::mass_matrix_type;
  }
  /*
    requirements on the nested aliases
  */
  && ::pressio::ops::is_known_data_type<typename T::state_type>::value
  && ::pressio::ops::is_known_data_type<typename T::right_hand_side_type>::value
  && ::pressio::ops::is_known_data_type<typename T::jacobian_type>::value
  && ::pressio::ops::is_known_data_type<typename T::mass_matrix_type>::value
  && all_have_traits_and_same_scalar<
    typename T::state_type,
    typename T::right_hand_side_type,
    typename T::jacobian_type,
    typename T::mass_matrix_type>::value
  && std::convertible_to<
    typename T::independent_variable_type, scalar_trait_t<typename T::state_type>>
  /*
    compund requirements on the "create" methods
  */
  && requires(const T & A)
  {
    { A.createState()         } -> std::same_as<typename T::state_type>;
    { A.createRightHandSide() } -> std::same_as<typename T::right_hand_side_type>;
    { A.createJacobian()      } -> std::same_as<typename T::jacobian_type>;
    { A.createMassMatrix()    } -> std::same_as<typename T::mass_matrix_type>;
  }
  /*
    compund requirements on "evaluation" method:
    intentionally not lumped with the above one for these reasons:
    1. it makes sense logically to split them, since this depends on the above
    2. helps the compiler with early failure detection
  */
  && requires(const T & A,
	      const typename T::state_type & state,
	      const typename T::independent_variable_type & evalValue,
	      typename T::right_hand_side_type & rhs,
	      typename T::jacobian_type & jac,
	      typename T::mass_matrix_type & M,
	      bool computeJacobian)
  {
    { A(state, evalValue, rhs, M, jac, computeJacobian) } -> std::same_as<void>;
  };

}} // end namespace pressio::ode




/* leave some white space on purpose so that
   if we make edits above, we don't have to change
   the line numbers included in the rst doc page */

#else


namespace pressio{ namespace ode{

template<class T, class enable = void>
struct SystemWithRhsJacobianMassMatrix : std::false_type{};

template<class T>
struct SystemWithRhsJacobianMassMatrix<
  T,
  mpl::enable_if_t<
    std::is_copy_constructible<T>::value
    //
    && ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_right_hand_side_typedef<T>::value
    && ::pressio::has_jacobian_typedef<T>::value
    && ::pressio::has_mass_matrix_typedef<T>::value
    //
    && ::pressio::ops::is_known_data_type<typename T::state_type>::value
    && ::pressio::ops::is_known_data_type<typename T::right_hand_side_type>::value
    && ::pressio::ops::is_known_data_type<typename T::jacobian_type>::value
    && ::pressio::ops::is_known_data_type<typename T::mass_matrix_type>::value
    && all_have_traits_and_same_scalar<
	 typename T::state_type,
	 typename T::right_hand_side_type,
	 typename T::jacobian_type,
	 typename T::mass_matrix_type>::value
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
    && ::pressio::ode::has_const_create_mass_matrix_method_return_result<
      T, typename T::mass_matrix_type >::value
    //
    && std::is_void<
      decltype(
	       std::declval<T const>()
	       (
		std::declval<typename T::state_type const&>(),
		std::declval<typename T::independent_variable_type const &>(),
		std::declval<typename T::right_hand_side_type &>(),
		std::declval<typename T::mass_matrix_type &>(),
		std::declval<typename T::jacobian_type &>(),
		std::declval<bool>()
	       )
	   )
      >::value
   >
  > : std::true_type{};

}} // end namespace pressio::ode

#endif //end PRESSIO_ENABLE_CXX20

#endif
