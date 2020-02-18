/*
//@HEADER
// ************************************************************************
//
// ode_is_legitimate_model_for_implicit_ode_arbitrary_stepper_with_user_defined_policies.hpp
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

#ifndef ODE_IS_LEGITIMATE_MODEL_FOR_IMPLICIT_ODE_ARBITRARY_STEPPER_WITH_USER_DEFINED_POLICIES_HPP_
#define ODE_IS_LEGITIMATE_MODEL_FOR_IMPLICIT_ODE_ARBITRARY_STEPPER_WITH_USER_DEFINED_POLICIES_HPP_

namespace pressio{ namespace ode{ namespace meta {

// when we have user-defined policies, a model does not necessarily have to have all the required
// time-discrete methods like we require for the arbitrary stepper with standard policies.
// This is because for custom policies we do not know what happens and it can very well
// be that the user handles things differently so the model can be anything.
// for instance, the user can define policies that call particular methods already existing in the model
// and therefore we do not care what the user does.
// Hence, when we have user-defined policies, we set the only requirement for the model to
// have the minimal set of typedefs

template<typename model_type, typename enable = void>
struct is_legitimate_model_for_implicit_ode_arbitrary_stepper_with_user_defined_policies
  : std::false_type{};

template<typename model_type>
struct is_legitimate_model_for_implicit_ode_arbitrary_stepper_with_user_defined_policies<
  model_type,
  mpl::enable_if_t<
    ::pressio::ode::meta::ode_model_has_all_needed_typedefs_for_implicit_ode_arbitrary_stepper<model_type>::value
    >
  > : std::true_type{};

}}} // namespace pressio::ode::meta
#endif
