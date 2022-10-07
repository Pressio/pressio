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

#ifndef ODE_STEPPERS_FULLY_DISCRETE_SYSTEM_HPP_
#define ODE_STEPPERS_FULLY_DISCRETE_SYSTEM_HPP_

#include "./predicates/ode_has_const_discrete_residual_jacobian_method.hpp"

#ifdef PRESSIO_ENABLE_CXX20


// this is here so that we can clearly show it in the doc via rst literal include directive
namespace pressio{ namespace ode{

template <class T, int NumStates>
concept FullyDiscreteSystemWithJacobian =
  std::copy_constructible<T>
  /* must have nested aliases */
  && requires(){
    typename T::independent_variable_type;
    typename T::state_type;
    typename T::discrete_residual_type;
    typename T::discrete_jacobian_type;
  }
  /*
    requirements on the nested aliases
  */
  && ::pressio::ops::is_known_data_type<typename T::state_type>::value
  && ::pressio::ops::is_known_data_type<typename T::discrete_residual_type>::value
  && ::pressio::ops::is_known_data_type<typename T::discrete_jacobian_type>::value
  && all_have_traits_and_same_scalar<
    typename T::state_type,
    typename T::discrete_residual_type,
    typename T::discrete_jacobian_type>::value
  && std::convertible_to<
    typename T::independent_variable_type, scalar_trait_t<typename T::state_type>>
  /*
    compund requirements on the "create" methods
  */
  && requires(const T & A)
  {
    { A.createState()            } -> std::same_as<typename T::state_type>;
    { A.createDiscreteResidual() } -> std::same_as<typename T::discrete_residual_type>;
    { A.createDiscreteJacobian() } -> std::same_as<typename T::discrete_jacobian_type>;
  }
  /*
    compund requirements on "evaluation" method:
    intentionally not lumped with the above one for these reasons:
    1. it makes sense logically to split them, since this depends on the above
    2. helps the compiler with early failure detection
  */
  && ::pressio::ode::has_const_discrete_residual_jacobian_method<
    T, NumStates,
    typename ::pressio::ode::StepCount::value_type,
    typename T::independent_variable_type,
    typename T::state_type,
    typename T::discrete_residual_type,
    typename T::discrete_jacobian_type>::value;

}} // end namespace pressio::ode





/* leave some white space on purpose so that
   if we make edits above, we don't have to change
   the line numbers included in the rst doc page */

#else

#include "./predicates/ode_has_const_create_state_method_return_result.hpp"
#include "./predicates/ode_has_const_create_discrete_residual_method_return_result.hpp"
#include "./predicates/ode_has_const_create_discrete_jacobian_method_return_result.hpp"

namespace pressio{ namespace ode{

template<class T, int NumStates,class enable = void>
struct FullyDiscreteSystemWithJacobian : std::false_type{};

template<class T, int NumStates>
struct FullyDiscreteSystemWithJacobian<
  T, NumStates,
  mpl::enable_if_t<
    std::is_copy_constructible<T>::value
    //
    && ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_discrete_residual_typedef<T>::value
    && ::pressio::has_discrete_jacobian_typedef<T>::value
    //
    && ::pressio::ops::is_known_data_type<typename T::state_type>::value
    && ::pressio::ops::is_known_data_type<typename T::discrete_residual_type>::value
    && ::pressio::ops::is_known_data_type<typename T::discrete_jacobian_type>::value
    && all_have_traits_and_same_scalar<
      typename T::state_type,
      typename T::discrete_residual_type,
      typename T::discrete_jacobian_type>::value
    && std::is_convertible<
      typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type>>::value
    //
    // creation methods
    && ::pressio::ode::has_const_create_state_method_return_result<
      T, typename T::state_type >::value
    && ::pressio::ode::has_const_create_discrete_residual_method_return_result<
      T, typename T::discrete_residual_type>::value
    && ::pressio::ode::has_const_create_discrete_jacobian_method_return_result<
      T, typename T::discrete_jacobian_type>::value
    //
    && ::pressio::ode::has_const_discrete_residual_jacobian_method<
      T, NumStates,
      typename ::pressio::ode::StepCount::value_type,
      typename T::independent_variable_type, typename T::state_type,
      typename T::discrete_residual_type, typename T::discrete_jacobian_type
      >::value
    >
  > : std::true_type{};

}} // end namespace pressio::ode
#endif //end PRESSIO_ENABLE_CXX20

#endif
