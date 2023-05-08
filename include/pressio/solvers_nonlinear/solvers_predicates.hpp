/*
//@HEADER
// ************************************************************************
//
// solvers_has_const_create_gradient_method_return_result.hpp
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

#ifndef SOLVERS_NONLINEAR_CONCEPTS_SOLVERS_PREDICATES_HPP_
#define SOLVERS_NONLINEAR_CONCEPTS_SOLVERS_PREDICATES_HPP_

#ifdef PRESSIO_ENABLE_CXX17
#include <optional>
#endif

namespace pressio{ namespace nonlinearsolvers{

template <class T, class StateType, class = void>
struct has_const_create_state_method_return_result
  : std::false_type{};

template <class T, class StateType>
struct has_const_create_state_method_return_result<
  T, StateType,
  ::pressio::mpl::enable_if_t<
    mpl::is_same<
      StateType,
      decltype(std::declval<T const>().createState())
      >::value
    >
  > : std::true_type{};

template<typename T, typename ResidualType, typename Enable = void>
struct has_const_create_residual_method_return_result : std::false_type{};

template<typename T, typename ResidualType>
struct has_const_create_residual_method_return_result
<T, ResidualType,
 ::pressio::mpl::enable_if_t<
   ::pressio::mpl::is_same<
     ResidualType,
     decltype( std::declval<T const>().createResidual() )
     >::value
   >
 > : std::true_type{};

template <
  class T,
  class StateType,
  class ResidualType,
  class = void
  >
struct has_const_residual_method_accept_state_result_return_void
  : std::false_type{};

template <
  class T,
  class StateType,
  class ResidualType
  >
struct has_const_residual_method_accept_state_result_return_void<
  T, StateType, ResidualType,
  mpl::enable_if_t<
    std::is_void<
      decltype(
         std::declval<T const>().residual
	 (
	    std::declval<StateType const &>(),
	    std::declval<ResidualType &>()
          )
         )
      >::value
    >
  > : std::true_type{};

template<class T, class JacobianType, class enable = void>
struct has_const_create_jacobian_method_return_result : std::false_type{};

template<class T, class JacobianType>
struct has_const_create_jacobian_method_return_result
<T, JacobianType,
 ::pressio::mpl::enable_if_t<
   ::pressio::mpl::is_same<
     JacobianType,
     decltype( std::declval<T const>().createJacobian() )
     >::value
   >
 > : std::true_type{};

template <
  class T,
  class StateType,
  class JacobianType,
  class = void
  >
struct has_const_jacobian_method_accept_state_result_return_void
  : std::false_type{};

template <
  class T,
  class StateType,
  class JacobianType
  >
struct has_const_jacobian_method_accept_state_result_return_void<
  T, StateType, JacobianType,
  mpl::enable_if_t<
    std::is_void<
      decltype(
         std::declval<T const>().jacobian
            (
              std::declval<StateType const &>(),
              std::declval<JacobianType &>()
            )
         )
      >::value
    >
  > : std::true_type{};


template <
  class T,
  class StateType,
  class ResidualType,
  class JacobianType,
  class = void
  >
struct has_const_residualandjacobian_method_accept_state_result_return_void
  : std::false_type{};

template <
  class T,
  class StateType,
  class ResidualType,
  class JacobianType
  >
struct has_const_residualandjacobian_method_accept_state_result_return_void<
  T, StateType, ResidualType, JacobianType,
  mpl::enable_if_t<
    std::is_void<
      decltype(
         std::declval<T const>().residualAndJacobian
            (
              std::declval<StateType const &>(),
              std::declval<ResidualType &>(),
#ifdef PRESSIO_ENABLE_CXX17
	      std::declval<std::optional<JacobianType*>>()
#else
	      std::declval<JacobianType*>()
#endif
            )
         )
      >::value
    >
  > : std::true_type{};

}} // namespace pressio::solvers
#endif  // SOLVERS_NONLINEAR_CONCEPTS_SOLVERS_PREDICATES_HPP_
