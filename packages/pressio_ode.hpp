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
#include "pressio_qr.hpp"
#include "pressio_svd.hpp"
#include "pressio_solvers.hpp"

#include "ode/src/ode_ConfigDefs.hpp"

#include "ode/src/ode_stepper_tags.hpp"
#include "ode/src/ode_types.hpp"
#include "ode/src/ode_fwd.hpp"

#include "ode/src/ode_aux_states_container.hpp"
#include "ode/src/ode_velocities_container.hpp"
#include "ode/src/ode_system_wrapper.hpp"

#include "ode/src/meta/ode_has_state_typedef.hpp"
#include "ode/src/meta/ode_has_velocity_typedef.hpp"
#include "ode/src/meta/ode_is_legitimate_guesser.hpp"
#include "ode/src/meta/ode_is_legitimate_time_step_size_setter.hpp"
#include "ode/src/meta/collector/ode_is_legitimate_collector.hpp"
#include "ode/src/meta/ode_model_has_all_needed_velocity_methods.hpp"

// --------------
// explicit
// --------------
#include "ode/src/explicit/meta/ode_is_legitimate_explicit_state_type.hpp"
#include "ode/src/explicit/meta/ode_is_legitimate_explicit_velocity_type.hpp"
#include "ode/src/explicit/meta/ode_is_legitimate_model_for_explicit_ode.hpp"
#include "ode/src/explicit/meta/ode_is_valid_user_defined_ops_for_explicit_euler.hpp"
#include "ode/src/explicit/meta/ode_is_valid_user_defined_ops_for_explicit_rk4.hpp"
#include "ode/src/explicit/meta/ode_is_valid_user_defined_ops_for_explicit_ode.hpp"

#include "ode/src/explicit/policies/ode_explicit_velocity_policy_base.hpp"
#include "ode/src/explicit/policies/ode_is_legitimate_explicit_velocity_policy.hpp"
#include "ode/src/explicit/policies/ode_explicit_velocity_standard_policy.hpp"
#include "ode/src/explicit/policies/ode_is_explicit_euler_velocity_standard_policy.hpp"
#include "ode/src/explicit/policies/ode_is_explicit_runge_kutta4_velocity_standard_policy.hpp"

#include "ode/src/explicit/steppers/ode_explicit_stepper_base.hpp"
#include "ode/src/explicit/steppers/ode_explicit_stepper.hpp"
#include "ode/src/explicit/steppers/ode_explicit_stepper_traits_euler.hpp"
#include "ode/src/explicit/steppers/ode_explicit_stepper_traits_rk4.hpp"

// --------------
// implicit
// --------------
#include "ode/src/implicit/ode_implicit_constants.hpp"
#include "ode/src/implicit/ode_time_discrete_residual.hpp"
#include "ode/src/implicit/ode_time_discrete_jacobian.hpp"

#include "ode/src/implicit/meta/valid_ud_ops/ode_is_valid_user_defined_ops_for_implicit_bdf2.hpp"
#include "ode/src/implicit/meta/valid_ud_ops/ode_is_valid_user_defined_ops_for_implicit_euler.hpp"
#include "ode/src/implicit/meta/valid_ud_ops/ode_is_valid_user_defined_ops_for_implicit_ode.hpp"

#include "ode/src/implicit/meta/ode_has_residual_typedef.hpp"
#include "ode/src/implicit/meta/ode_has_jacobian_typedef.hpp"
#include "ode/src/implicit/meta/ode_is_stepper_order_setter.hpp"
#include "ode/src/implicit/meta/ode_is_stepper_total_n_states_setter.hpp"
#include "ode/src/implicit/meta/ode_is_legitimate_implicit_state_type.hpp"
#include "ode/src/implicit/meta/ode_is_legitimate_implicit_residual_type.hpp"
#include "ode/src/implicit/meta/ode_is_legitimate_implicit_jacobian_type.hpp"
#include "ode/src/implicit/meta/ode_is_legitimate_solver_for_implicit_stepper.hpp"
#include "ode/src/implicit/meta/ode_is_legitimate_auxiliary_stepper_for_implicit_ode.hpp"
#include "ode/src/implicit/meta/ode_implicit_stepper_stencil_needs_previous_states_and_velocities.hpp"
#include "ode/src/implicit/meta/ode_implicit_stepper_stencil_needs_previous_states_only.hpp"

// metaf for model for arbitrary stepper
#include "ode/src/implicit/meta/admissible_model_arbitrary_stepper/ode_model_has_all_needed_typedefs_for_implicit_ode_arbitrary_stepper.hpp"
#include "ode/src/implicit/meta/admissible_model_arbitrary_stepper/ode_has_time_discrete_jacobian_method_accepting_n_states_returning_void.hpp"
#include "ode/src/implicit/meta/admissible_model_arbitrary_stepper/ode_has_create_time_discrete_jacobian_object_method_returning_non_void.hpp"
#include "ode/src/implicit/meta/admissible_model_arbitrary_stepper/ode_has_needed_time_discrete_jacobian_methods.hpp"
#include "ode/src/implicit/meta/admissible_model_arbitrary_stepper/ode_has_time_discrete_residual_method_accepting_n_states_returning_void.hpp"
#include "ode/src/implicit/meta/admissible_model_arbitrary_stepper/ode_has_create_time_discrete_residual_object_method_returning_non_void.hpp"
#include "ode/src/implicit/meta/admissible_model_arbitrary_stepper/ode_has_needed_time_discrete_residual_methods.hpp"
#include "ode/src/implicit/meta/admissible_model_arbitrary_stepper/ode_is_legitimate_model_for_implicit_ode_arbitrary_stepper_with_standard_policies.hpp"
#include "ode/src/implicit/meta/admissible_model_arbitrary_stepper/ode_is_legitimate_model_for_implicit_ode_arbitrary_stepper_with_standard_res_ud_jac_policies.hpp"
#include "ode/src/implicit/meta/admissible_model_arbitrary_stepper/ode_is_legitimate_model_for_implicit_ode_arbitrary_stepper_with_ud_res_standard_jac_policies.hpp"
#include "ode/src/implicit/meta/admissible_model_arbitrary_stepper/ode_is_legitimate_model_for_implicit_ode_arbitrary_stepper_with_user_defined_policies.hpp"
#include "ode/src/implicit/meta/admissible_model_arbitrary_stepper/ode_model_is_compatible_with_policies_types_for_implicit_ode_arbitrary_stepper.hpp"
#include "ode/src/implicit/meta/admissible_model_arbitrary_stepper/ode_is_legitimate_model_for_implicit_ode_arbitrary_stepper.hpp"

// metaf for model for regular stepper
#include "ode/src/implicit/meta/admissible_model_regular_stepper/ode_model_has_all_needed_typedefs_for_implicit_ode_regular_stepper.hpp"
#include "ode/src/implicit/meta/admissible_model_regular_stepper/ode_has_jacobian_method_callable_with_two_args.hpp"
#include "ode/src/implicit/meta/admissible_model_regular_stepper/ode_has_jacobian_method_callable_with_three_args.hpp"
#include "ode/src/implicit/meta/admissible_model_regular_stepper/ode_model_has_all_needed_jacobian_methods.hpp"
#include "ode/src/implicit/meta/admissible_model_regular_stepper/ode_is_legitimate_model_for_implicit_ode_regular_stepper_with_standard_policies.hpp"
#include "ode/src/implicit/meta/admissible_model_regular_stepper/ode_is_legitimate_model_for_implicit_ode_regular_stepper_with_user_defined_policies.hpp"
#include "ode/src/implicit/meta/admissible_model_regular_stepper/ode_is_legitimate_model_for_implicit_ode_regular_stepper_with_standard_res_ud_jac_policies.hpp"
#include "ode/src/implicit/meta/admissible_model_regular_stepper/ode_is_legitimate_model_for_implicit_ode_regular_stepper_with_ud_res_standard_jac_policies.hpp"
#include "ode/src/implicit/meta/admissible_model_regular_stepper/ode_model_is_compatible_with_policies_types_for_implicit_ode_regular_stepper.hpp"
#include "ode/src/implicit/meta/admissible_model_regular_stepper/ode_is_legitimate_model_for_implicit_ode_regular_stepper.hpp"

#include "ode/src/implicit/meta/ode_is_legitimate_model_for_implicit_ode.hpp"

// implicit policies
#include "ode/src/implicit/policies/standard/ode_implicit_residual_standard_policy.hpp"
#include "ode/src/implicit/policies/standard/ode_implicit_residual_standard_policy_for_arbitrary_stepper.hpp"
#include "ode/src/implicit/policies/standard/ode_implicit_residual_standard_policy_pybind11.hpp"
#include "ode/src/implicit/policies/standard/ode_implicit_jacobian_standard_policy.hpp"
#include "ode/src/implicit/policies/standard/ode_implicit_jacobian_standard_policy_for_arbitrary_stepper.hpp"
#include "ode/src/implicit/policies/standard/ode_implicit_jacobian_standard_policy_pybind11.hpp"

#include "ode/src/implicit/policies/meta/ode_is_implicit_residual_standard_policy.hpp"
#include "ode/src/implicit/policies/meta/ode_is_implicit_jacobian_standard_policy.hpp"
#include "ode/src/implicit/policies/meta/ode_is_legitimate_implicit_residual_policy.hpp"
#include "ode/src/implicit/policies/meta/ode_is_legitimate_implicit_jacobian_policy.hpp"
#include "ode/src/implicit/policies/meta/ode_find_if_legitimate_implicit_jacobian_policy.hpp"
#include "ode/src/implicit/policies/meta/ode_find_if_legitimate_implicit_residual_policy.hpp"
#include "ode/src/implicit/policies/meta/ode_is_legitimate_residual_policy_for_implicit_arbitrary_stepper.hpp"
#include "ode/src/implicit/policies/meta/ode_is_legitimate_jacobian_policy_for_implicit_arbitrary_stepper.hpp"
#include "ode/src/implicit/policies/meta/ode_find_if_legitimate_jacobian_policy_for_implicit_arbitrary_stepper.hpp"
#include "ode/src/implicit/policies/meta/ode_find_if_legitimate_residual_policy_for_implicit_arbitrary_stepper.hpp"

// steppers
#include "ode/src/implicit/steppers/ode_implicit_stepper_base.hpp"
#include "ode/src/implicit/steppers/ode_implicit_stepper_arbitrary.hpp"
#include "ode/src/implicit/steppers/ode_implicit_stepper_bdf2.hpp"
#include "ode/src/implicit/steppers/ode_implicit_stepper_euler.hpp"
#include "ode/src/implicit/steppers/ode_implicit_stepper_traits_helpers.hpp"
#include "ode/src/implicit/steppers/ode_implicit_stepper_traits_static_checks.hpp"
#include "ode/src/implicit/steppers/ode_implicit_stepper_traits_arbitrary.hpp"
#include "ode/src/implicit/steppers/ode_implicit_stepper_traits_bdf2.hpp"
#include "ode/src/implicit/steppers/ode_implicit_stepper_traits_euler.hpp"

// --------------
// integrators
// --------------
#include "ode/src/integrators/ode_integrate_n_steps_explicit.hpp"
#include "ode/src/integrators/ode_integrate_n_steps_implicit_arbitrary_step_size.hpp"
#include "ode/src/integrators/ode_integrate_n_steps_implicit_constant_step_size.hpp"
#include "ode/src/integrators/ode_integrate_to_target_time_implicit_arbitrary_step_size.hpp"

#endif
