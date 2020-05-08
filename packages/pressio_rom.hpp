/*
//@HEADER
// ************************************************************************
//
// pressio_rom.hpp
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

#ifndef PRESSIO_ROM_HPP_
#define PRESSIO_ROM_HPP_

#include "pressio_mpl.hpp"
#include "pressio_utils.hpp"
#include "pressio_containers.hpp"
#include "pressio_ops.hpp"
#include "pressio_qr.hpp"
#include "pressio_svd.hpp"
#include "pressio_optimizers.hpp"
#include "pressio_solvers.hpp"
#include "pressio_ode.hpp"

#include "rom/src/rom_ConfigDefs.hpp"

#include "rom/src/rom_static_container_fom_states.hpp"
#include "rom/src/utils/rom_utils_set_gen_coordinates.hpp"
#include "rom/src/meta/rom_has_dense_matrix_typedef.hpp"

//----------------
// decorators
//----------------
#include "rom/src/decorators/rom_preconditioned_decorator_residual.hpp"
#include "rom/src/decorators/rom_preconditioned_decorator_jacobian.hpp"
#include "rom/src/decorators/rom_mask_decorator_residual.hpp"
#include "rom/src/decorators/rom_mask_decorator_jacobian.hpp"

//----------------
// decoder classes
//----------------
#include "rom/src/meta/decoder/rom_is_legitimate_custom_ops_for_linear_decoder.hpp"
#include "rom/src/meta/decoder/rom_is_legitimate_decoder_type.hpp"
#include "rom/src/decoder/rom_linear_decoder.hpp"

//----------------------
// fom-querying policies
//----------------------
#include "rom/src/fom_querying_policies/rom_query_fom_velocity_unsteady.hpp"
#include "rom/src/fom_querying_policies/rom_query_fom_velocity_steady.hpp"
#include "rom/src/fom_querying_policies/rom_query_fom_apply_jacobian_unsteady.hpp"
#include "rom/src/fom_querying_policies/rom_query_fom_apply_jacobian_steady.hpp"
#include "rom/src/fom_querying_policies/rom_query_fom_apply_time_discrete_jacobian.hpp"
#include "rom/src/fom_querying_policies/rom_query_fom_time_discrete_residual.hpp"

//----------------
// classes for fom state reconstructor
//----------------
#include "rom/src/meta/rom_is_legitimate_custom_ops_for_fom_state_reconstructor.hpp"
#include "rom/src/fom_state_reconstructor/rom_reconstructor_fom_state.hpp"

//----------------
// steady LSPG
//----------------
#include "rom/src/lspg_steady/rom_lspg_steady_residual_policy.hpp"
#include "rom/src/lspg_steady/rom_lspg_steady_jacobian_policy.hpp"
#include "rom/src/lspg_steady/rom_lspg_steady_system.hpp"
#include "rom/src/lspg_steady/problem_traits/rom_lspg_steady_common_traits.hpp"
#include "rom/src/lspg_steady/problem_traits/rom_lspg_steady_default_problem_traits.hpp"
#include "rom/src/lspg_steady/problem_traits/rom_lspg_steady_preconditioned_problem_traits.hpp"
#include "rom/src/lspg_steady/rom_lspg_steady_problem_generator.hpp"
#include "rom/src/lspg_steady/rom_lspg_steady_api_aliases.hpp"

//----------------
// unsteady LSPG
//----------------
// metaf for velocity api
#include "rom/src/meta/lspg_velocity_api/rom_is_legitimate_custom_ops_for_unsteady_lspg_velocity_api.hpp"
#include "rom/src/meta/lspg_velocity_api/rom_has_apply_jacobian_method_callable_with_three_args_for_unsteady.hpp"
#include "rom/src/meta/lspg_velocity_api/rom_has_apply_jacobian_method_callable_with_four_args_for_unsteady.hpp"
#include "rom/src/meta/lspg_velocity_api/rom_model_has_needed_apply_jacobian_methods_for_unsteady.hpp"
#include "rom/src/meta/lspg_velocity_api/rom_model_has_needed_velocity_methods.hpp"
#include "rom/src/meta/lspg_velocity_api/rom_model_meets_velocity_api_for_unsteady_lspg.hpp"

// metaf for residual api
#include "rom/src/meta/lspg_residual_api/rom_is_legitimate_custom_ops_for_unsteady_lspg_residual_api.hpp"
#include "rom/src/meta/lspg_residual_api/rom_has_apply_time_discrete_jacobian_method_accepting_n_states_returning_void.hpp"
#include "rom/src/meta/lspg_residual_api/rom_has_create_apply_time_discrete_jacobian_object_method_returning_non_void.hpp"
#include "rom/src/meta/lspg_residual_api/rom_model_has_needed_apply_time_discrete_jacobian_methods.hpp"
#include "rom/src/meta/lspg_residual_api/rom_model_has_needed_typedefs_for_unsteady_lspg_residual_api.hpp"
#include "rom/src/meta/lspg_residual_api/rom_model_meets_residual_api_for_unsteady_lspg.hpp"

// metaf for checking valid model for unsteady lspg
#include "rom/src/meta/rom_is_legitimate_model_for_unsteady_lspg.hpp"
// lspg problems
#include "rom/src/lspg_unsteady/rom_lspg_unsteady_problem_default.hpp"
#include "rom/src/lspg_unsteady/rom_lspg_unsteady_problem_masked.hpp"
#include "rom/src/lspg_unsteady/rom_lspg_unsteady_problem_preconditioned.hpp"
#include "rom/src/lspg_unsteady/rom_lspg_unsteady_problem_generator.hpp"

//----------------
// galerkin
//----------------
// metaf for checking ops
#include "rom/src/meta/galerkin_velocity_api/rom_is_legitimate_custom_ops_for_galerkin_velocity_api.hpp"
// meta for checking api
#include "rom/src/meta/galerkin_velocity_api/rom_model_meets_velocity_api_for_galerkin.hpp"
#include "rom/src/meta/galerkin_residual_api/rom_model_meets_residual_api_for_galerkin.hpp"
#include "rom/src/meta/rom_is_legitimate_model_for_galerkin.hpp"
// problem classes
#include "rom/src/galerkin/rom_galerkin_problem_default.hpp"
#include "rom/src/galerkin/rom_galerkin_problem_generator.hpp"

//----------------
// wls
//----------------
#include "rom/src/meta/wls_velocity_api/rom_model_meets_velocity_api_for_wls.hpp"
#include "rom/src/meta/wls_residual_api/rom_model_meets_residual_api_for_wls.hpp"
#include "rom/src/wls/rom_wls_types.hpp"
#include "rom/src/wls/rom_wls_preconditioners.hpp"
#include "rom/src/wls/meta/rom_wls_is_legitimate_preconditioner_type.hpp"
#include "rom/src/wls/rom_wls_jacobian_updating_tag.hpp"
#include "rom/src/wls/meta/rom_wls_is_legitimate_jacobian_updating_tag.hpp"
#include "rom/src/wls/rom_wls_jacobians_container.hpp"
#include "rom/src/wls/time_schemes/rom_wls_implicit_euler.hpp"
#include "rom/src/wls/time_schemes/rom_wls_bdf2.hpp"
#include "rom/src/wls/time_schemes/rom_wls_select_timescheme_helper.hpp"
#include "rom/src/wls/policies/rom_wls_hessian_and_gradient_sequential_policy.hpp"
#include "rom/src/wls/apis/rom_wls_hessian_gradient_system_api.hpp"

//----------------
// Common tags
//----------------

#endif
