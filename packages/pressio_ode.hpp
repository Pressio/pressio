/*
//@HEADER
// ************************************************************************
//
// pressio_ode.hpp
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

#ifndef PRESSIO_ODE_HPP_
#define PRESSIO_ODE_HPP_

#include "pressio_mpl.hpp"
#include "pressio_utils.hpp"
#include "pressio_containers.hpp"
#include "pressio_ops.hpp"
#include "pressio_qr.hpp"
#include "pressio_svd.hpp"
#include "pressio_solvers.hpp"

#include "ode/src/ode_types.hpp"
#include "ode/src/ode_fwd.hpp"
#include "ode/src/ode_aux_states_container.hpp"
#include "ode/src/ode_velocities_container.hpp"
#include "ode/src/ode_system_wrapper.hpp"

#include "ode/src/predicates/typedefs/ode_has_state_typedef.hpp"
#include "ode/src/predicates/typedefs/ode_has_velocity_typedef.hpp"
#include "ode/src/predicates/typedefs/ode_has_residual_typedef.hpp"
#include "ode/src/predicates/typedefs/ode_has_jacobian_typedef.hpp"

#include "ode/src/will_be_concepts/collector/ode_collector.hpp"
#include "ode/src/will_be_concepts/ode_guesser.hpp"
#include "ode/src/will_be_concepts/ode_legitimate_time_step_size_setter.hpp"

// --------------
// explicit
// --------------
#include "ode/src/explicit/ode_explicit_stepper_tags.hpp"

#include "ode/src/will_be_concepts/ode_explicit_state_type.hpp"
#include "ode/src/will_be_concepts/ode_explicit_velocity_type.hpp"
#include "ode/src/predicates/velocity_methods/ode_has_const_create_velocity_method_return_result.hpp"
#include "ode/src/predicates/velocity_methods/ode_has_const_velocity_method_accept_state_time_result_return_void.hpp"
#include "ode/src/will_be_concepts/policies/ode_admissible_explicit_velocity_policy.hpp"
#include "ode/src/will_be_concepts/system/ode_admissible_system_explicit_ode.hpp"
#include "ode/src/will_be_concepts/user_defined_ops/ode_admissible_user_defined_ops_for_explicit_euler.hpp"
#include "ode/src/will_be_concepts/user_defined_ops/ode_admissible_user_defined_ops_for_explicit_rk4.hpp"
#include "ode/src/will_be_concepts/system/ode_admissible_system_explicit_ode.hpp"
#include "ode/src/will_be_concepts/policies/ode_admissible_explicit_velocity_policy.hpp"

#include "ode/src/explicit/ode_explicit_velocity_standard_policy.hpp"
#include "ode/src/explicit/ode_explicit_stepper_base.hpp"
#include "ode/src/explicit/ode_explicit_stepper.hpp"

// --------------
// implicit
// --------------
#include "ode/src/implicit/ode_implicit_stepper_tags.hpp"

#include "ode/src/predicates/ode_is_stepper_total_n_states_setter.hpp"
#include "ode/src/predicates/ode_is_stepper_order_setter.hpp"
#include "ode/src/predicates/time_discrete_residual_methods/ode_has_const_create_td_residual_method_return_result.hpp"
#include "ode/src/predicates/time_discrete_residual_methods/ode_has_const_td_residual_method_accept_step_time_dt_result_norm_states_return_void.hpp"
#include "ode/src/predicates/time_discrete_jacobian_methods/ode_has_const_create_td_jacobian_method_return_result.hpp"
#include "ode/src/predicates/time_discrete_jacobian_methods/ode_has_const_time_discrete_jacobian_method_accepting_n_states_returning_void.hpp"

#include "ode/src/will_be_concepts/ode_implicit_state_type.hpp"
#include "ode/src/will_be_concepts/ode_implicit_residual_type.hpp"
#include "ode/src/will_be_concepts/ode_implicit_jacobian_type.hpp"
#include "ode/src/predicates/jacobian_methods/ode_has_const_create_jacobian_method_return_result.hpp"
#include "ode/src/predicates/jacobian_methods/ode_has_const_jacobian_method_accept_state_time_result_return_void.hpp"
#include "ode/src/will_be_concepts/ode_legitimate_solver_for_implicit_stepper.hpp"
#include "ode/src/will_be_concepts/ode_admissible_auxiliary_stepper_for_implicit_ode.hpp"
#include "ode/src/will_be_concepts/policies/ode_admissible_implicit_residual_policy_regular_stepper.hpp"
#include "ode/src/will_be_concepts/policies/ode_admissible_implicit_jacobian_policy_regular_stepper.hpp"
#include "ode/src/will_be_concepts/policies/ode_admissible_implicit_residual_policy_arbitrary_stepper.hpp"
#include "ode/src/will_be_concepts/policies/ode_admissible_implicit_jacobian_policy_arbitrary_stepper.hpp"
#include "ode/src/will_be_concepts/system/ode_admissible_system_implicit_ode_regular_stepper.hpp"
#include "ode/src/will_be_concepts/system/ode_admissible_system_implicit_ode_arbitrary_stepper.hpp"
#include "ode/src/will_be_concepts/system/ode_admissible_system_implicit_ode.hpp"

#include "ode/src/implicit/ode_implicit_constants.hpp"
#include "ode/src/implicit/impl/ode_time_discrete_residual_impl.hpp"
#include "ode/src/implicit/impl/ode_time_discrete_jacobian_impl.hpp"
#include "ode/src/implicit/standard_policies/ode_implicit_residual_standard_policy.hpp"
#include "ode/src/implicit/standard_policies/ode_implicit_residual_standard_policy_for_arbitrary_stepper.hpp"
#include "ode/src/implicit/standard_policies/ode_implicit_jacobian_standard_policy.hpp"
#include "ode/src/implicit/standard_policies/ode_implicit_jacobian_standard_policy_for_arbitrary_stepper.hpp"
#include "ode/src/implicit/ode_implicit_stepper_base.hpp"
#include "ode/src/implicit/ode_implicit_stepper.hpp"

// --------------
// integrators
// --------------
#include "ode/src/integrators/ode_integrate_n_steps_explicit.hpp"
#include "ode/src/integrators/ode_integrate_n_steps_implicit_arbitrary_step_size.hpp"
#include "ode/src/integrators/ode_integrate_n_steps_implicit_constant_step_size.hpp"
#include "ode/src/integrators/ode_integrate_to_target_time_implicit_arbitrary_step_size.hpp"

#endif
