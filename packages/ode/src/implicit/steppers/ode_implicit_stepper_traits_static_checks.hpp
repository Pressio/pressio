/*
//@HEADER
// ************************************************************************
//
// ode_implicit_stepper_traits_static_checks.hpp
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

#ifndef ODE_STEPPERS_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_TRAITS_STATIC_CHECKS_HPP_
#define ODE_STEPPERS_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_TRAITS_STATIC_CHECKS_HPP_

#include "../../ode_fwd.hpp"
#include "../meta/ode_is_valid_user_defined_ops_for_implicit_ode.hpp"
#include "../meta/ode_is_legitimate_model_for_implicit_ode.hpp"
#include "../meta/ode_is_legitimate_model_for_implicit_ode_arbitrary_stepper.hpp"
#include "../meta/ode_is_stepper_order_setter.hpp"
#include "../meta/ode_is_stepper_total_n_states_setter.hpp"

namespace pressio{ namespace ode{ namespace details{ namespace impl{

template <typename T, bool isStdPolicy>
struct CheckModelTimeDiscreteJacobianMethods;

template <typename T>
struct CheckModelTimeDiscreteJacobianMethods<T, false> {};

template <typename T>
struct CheckModelTimeDiscreteJacobianMethods<T, true>
{
  static_assert(::pressio::ode::meta::has_needed_time_discrete_jacobian_methods<
		T, types::step_t,
		typename T::scalar_type, typename T::state_type, typename T::jacobian_type
		>::value,
		"\nThe model type you passed to the Arbitrary implicit stepper \n \
does not have valid jacobian methods, see api for reference.");
};


}}}}//end namespace pressio::ode::details::impl
#endif