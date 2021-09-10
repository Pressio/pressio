/*
//@HEADER
// ************************************************************************
//
// ode_implicit_jacobian_policy.hpp
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

#ifndef ODE_IMPLICIT_CONSTRAINTS_ODE_IMPLICIT_JACOBIAN_POLICY_HPP_
#define ODE_IMPLICIT_CONSTRAINTS_ODE_IMPLICIT_JACOBIAN_POLICY_HPP_

namespace pressio{ namespace ode{

template<
  class T,
  class StateType,
  class ScalarType,
  class = void
  >
struct implicit_jacobian_policy : std::false_type{};

template<
  class T,
  class StateType,
  class ScalarType
  >
struct implicit_jacobian_policy<
  T, StateType, ScalarType,
  ::pressio::mpl::enable_if_t<
    ::pressio::has_jacobian_typedef<T>::value
    and
    std::is_same<
      typename T::jacobian_type,
      decltype(std::declval<T const>().create())
      >::value
    and
    //
    std::is_void<
      decltype
      (
       std::declval<T const>()
       (
	std::declval<StepScheme const &>(),
	std::declval<StateType const &>(),
	std::declval<ImplicitStencilStatesContainerDyn<StateType> const & >(),
	std::declval<ScalarType const &>(),
	std::declval<ScalarType const &>(),
	std::declval<int const &>(),
	std::declval<typename T::jacobian_type &>()
	)
       )
      >::value
    >
  > : std::true_type{};
//------------------------------------------------------------------

template<class T, class ... args>
using implicit_euler_jacobian_policy =
  implicit_jacobian_policy<T, args...>;

template<class T, class ... args>
using implicit_bdf2_jacobian_policy =
  implicit_jacobian_policy<T, args...>;

template<class T, class ... args>
using implicit_cranknicolson_jacobian_policy =
  implicit_jacobian_policy<T, args...>;

}} // namespace pressio::ode::constraints
#endif  // ODE_IMPLICIT_CONSTRAINTS_ODE_IMPLICIT_JACOBIAN_POLICY_HPP_
