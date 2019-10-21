/*
//@HEADER
// ************************************************************************
//
// ode_model_is_compatible_with_policies_for_implicit_ode.hpp
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

#ifndef ODE_MODEL_IS_COMPATIBLE_WITH_POLICIES_FOR_IMPLICIT_ODE_HPP_
#define ODE_MODEL_IS_COMPATIBLE_WITH_POLICIES_FOR_IMPLICIT_ODE_HPP_

#include "ode_is_legitimate_model_for_implicit_ode_regular_stepper_with_standard_policies.hpp"
#include "ode_is_legitimate_model_for_implicit_ode_regular_stepper_with_user_defined_policies.hpp"

namespace pressio{ namespace ode{ namespace meta {

template<
  typename model_type,
  bool residual_policy_is_standard,
  bool jacobian_policy_is_standard,
  typename enable = void
  >
struct model_is_compatible_with_policies_for_implicit_ode : std::false_type{};

// specialize when policies are standard
template<typename model_type>
struct model_is_compatible_with_policies_for_implicit_ode<
  model_type, true, true,
  mpl::enable_if_t<
    is_legitimate_model_for_implicit_ode_regular_stepper_with_standard_policies<model_type>::value
    >
  > : std::true_type{};

// specialize when policies are custom
template<typename model_type>
struct model_is_compatible_with_policies_for_implicit_ode<
  model_type, false, false,
  mpl::enable_if_t<
    is_legitimate_model_for_implicit_ode_regular_stepper_with_user_defined_policies<model_type>::value
    >
  > : std::true_type{};


}}} // namespace pressio::ode::meta
#endif
