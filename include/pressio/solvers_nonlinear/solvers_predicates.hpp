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

#ifndef SOLVERS_NONLINEAR_PREDICATES_HPP_
#define SOLVERS_NONLINEAR_PREDICATES_HPP_

namespace pressio{ namespace nonlinearsolvers{

template<typename T, typename GradientType, typename enable = void>
struct has_const_create_gradient_method_return_result : std::false_type{};

template<typename T, typename GradientType>
struct has_const_create_gradient_method_return_result
<T, GradientType,
 ::pressio::mpl::enable_if_t<
   ::pressio::mpl::is_same<
     GradientType,
     decltype( std::declval<T const>().createGradient() )
     >::value
   >
 > : std::true_type{};


template<typename T, typename HessianType, typename enable = void>
struct has_const_create_hessian_method_return_result : std::false_type{};

template<typename T, typename HessianType>
struct has_const_create_hessian_method_return_result
<T, HessianType,
 ::pressio::mpl::enable_if_t<
   ::pressio::mpl::is_same<
     HessianType,
     decltype( std::declval<T const>().createHessian() )
     >::value
   >
 > : std::true_type{};


template<typename T, typename JacobianType, typename enable = void>
struct has_const_create_jacobian_method_return_result : std::false_type{};

template<typename T, typename JacobianType>
struct has_const_create_jacobian_method_return_result
<T, JacobianType,
 ::pressio::mpl::enable_if_t<
   ::pressio::mpl::is_same<
     JacobianType,
     decltype( std::declval<T const>().createJacobian() )
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


template <
  typename T,
  typename StateType,
  typename GradientType,
  typename NormType,
  typename = void
  >
struct has_const_gradient_method_accept_state_result_norm_return_void
  : std::false_type{};

template <
  typename T,
  typename StateType,
  typename GradientType,
  typename NormType
  >
struct has_const_gradient_method_accept_state_result_norm_return_void<
  T, StateType, GradientType, NormType,
  mpl::enable_if_t<
    std::is_void<
      decltype(
         std::declval<T const>().gradient(
            std::declval<StateType const &>(),
            std::declval<GradientType &>(),
            ::pressio::Norm::Undefined,
            std::declval<NormType &>(),
      std::declval<bool>()
            )
         )
      >::value
    >
  > : std::true_type{};


template <
  typename T,
  typename StateType,
  typename HessianType,
  typename = void
  >
struct has_const_hessian_method_accept_state_result_return_void
  : std::false_type{};

template <
  typename T,
  typename StateType,
  typename HessianType
  >
struct has_const_hessian_method_accept_state_result_return_void<
  T, StateType, HessianType, 
  mpl::enable_if_t<
    std::is_void<
      decltype(
         std::declval<T const>().hessian(
            std::declval<StateType const &>(),
            std::declval<HessianType &>()
            )
         )
      >::value
    >
  > : std::true_type{};


template <
  typename T,
  typename StateType,
  typename HessianType,
  typename GradientType,
  typename NormType,
  typename = void
  >
struct has_const_hessianandgradient_method_accept_state_result_norm_return_void
  : std::false_type{};

template <
  typename T,
  typename StateType,
  typename HessianType,
  typename GradientType,
  typename NormType
  >
struct has_const_hessianandgradient_method_accept_state_result_norm_return_void<
  T, StateType, HessianType, GradientType, NormType,
  mpl::enable_if_t<
    std::is_void<
      decltype(
         std::declval<T const>().hessianAndGradient
          (
            std::declval<StateType const &>(),
            std::declval<HessianType &>(),
            std::declval<GradientType &>(),
            /* does not matter here what we pass, just to test */
            ::pressio::Norm::Undefined,
            std::declval<NormType &>(),
      std::declval<bool>()
          )
         )
      >::value
    >
  > : std::true_type{};


template <
  typename T,
  typename StateType,
  typename JacobianType,
  typename = void
  >
struct has_const_jacobian_method_accept_state_result_return_void
  : std::false_type{};

template <
  typename T,
  typename StateType,
  typename JacobianType
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
  typename T,
  typename StateType,
  typename ResidualType,
  typename = void
  >
struct has_const_residual_method_accept_state_result_return_void
  : std::false_type{};

template <
  typename T,
  typename StateType,
  typename ResidualType
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


template <
  typename T,
  typename StateType,
  typename ResidualType,
  typename JacobianType,
  typename = void
  >
struct has_const_residualandjacobian_method_accept_state_result_return_void
  : std::false_type{};

template <
  typename T,
  typename StateType,
  typename ResidualType,
  typename JacobianType
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
              std::declval<JacobianType &>(),
              std::declval<bool>() //bool for updating or not jacobian
            )
         )
      >::value
    >
  > : std::true_type{};


template <
  typename T,
  typename StateType,
  typename NormType,
  typename = void
  >
struct has_const_residualnorm_method_accept_state_norm_return_void
  : std::false_type{};

template <
  typename T,
  typename StateType,
  typename NormType
  >
struct has_const_residualnorm_method_accept_state_norm_return_void<
  T, StateType, NormType, 
  mpl::enable_if_t<
    std::is_void<
      decltype(
         std::declval<T const>().residualNorm
            (
              std::declval<StateType const &>(),
              ::pressio::Norm::Undefined,
              std::declval<NormType &>()
            )
         )
      >::value
    >
  > : std::true_type{};

    
}} // namespace pressio::solvers
#endif
