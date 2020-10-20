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

#include "rom/src/rom_fwd.hpp"

#include "rom/src/rom_manager_fom_states_static.hpp"
#include "rom/src/utils/rom_utils_set_gen_coordinates.hpp"
#include "rom/src/predicates/typedefs/rom_has_dense_matrix_typedef.hpp"
#include "rom/src/will_be_concepts/rom_fom_state.hpp"
#include "rom/src/will_be_concepts/rom_rom_state.hpp"

// custom ops
#include "rom/src/will_be_concepts/custom_ops/rom_custom_ops_for_linear_decoder.hpp"
#include "rom/src/will_be_concepts/custom_ops/rom_custom_ops_for_fom_state_reconstructor.hpp"
#include "rom/src/will_be_concepts/custom_ops/rom_custom_ops_galerkin_continuous_time.hpp"
#include "rom/src/will_be_concepts/custom_ops/rom_custom_ops_lspg_continuous_time.hpp"
#include "rom/src/will_be_concepts/custom_ops/rom_custom_ops_lspg_discrete_time.hpp"

// decoder
#include "rom/src/predicates/decoder/rom_has_const_apply_mapping_accept_operand_result_return_void.hpp"
#include "rom/src/predicates/decoder/rom_has_const_get_reference_to_jacobian.hpp"
#include "rom/src/predicates/decoder/rom_has_const_update_jacobian_method_accept_operand_return_void.hpp"
#include "rom/src/will_be_concepts/decoder/rom_admissible_decoder.hpp"
#include "rom/src/will_be_concepts/decoder/rom_decoder_jacobian.hpp"
#include "rom/src/decoder/rom_linear_decoder.hpp"
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
#include "rom/src/decoder/rom_py_decoder.hpp"
#endif

// fom state reconstructor
#include "rom/src/fom_state_reconstructor/rom_reconstructor_fom_state.hpp"

// decorators
#include "rom/src/decorators/rom_preconditioned_decorator_residual.hpp"
#include "rom/src/decorators/rom_preconditioned_decorator_jacobian.hpp"
#include "rom/src/decorators/rom_masked_residual_policy.hpp"
#include "rom/src/decorators/rom_masked_jacobian_policy.hpp"

//---------------------------------
// system predicates and concepts
//---------------------------------
#include "rom/src/predicates/masking_methods/rom_has_const_apply_mask_method_accept_operand_result_return_void.hpp"
#include "rom/src/predicates/masking_methods/rom_has_const_apply_mask_method_accept_operand_time_result_return_void.hpp"
#include "rom/src/predicates/masking_methods/rom_has_const_create_apply_mask_result_method_accept_operand_return_result.hpp"
#include "rom/src/predicates/apply_discrete_time_jacobian_methods/rom_has_const_apply_discrete_time_jacobian_method_accept_step_time_dt_operand_result_n_states_returning_void.hpp"
#include "rom/src/predicates/apply_discrete_time_jacobian_methods/rom_has_const_create_apply_discrete_time_jacobian_result_method_accept_operand_return_result.hpp"
#include "rom/src/predicates/apply_jacobian_methods/rom_has_const_create_apply_jacobian_result_method_accept_operand_return_result.hpp"
#include "rom/src/predicates/apply_jacobian_methods/rom_has_const_apply_jacobian_method_accept_state_operand_time_result_return_void.hpp"
#include "rom/src/predicates/preconditioning_methods/rom_has_const_apply_preconditioner_method_accept_state_time_result_return_void.hpp"
#include "rom/src/predicates/residual_methods/rom_has_const_create_residual_method_return_result.hpp"
#include "rom/src/predicates/residual_methods/rom_has_const_residual_method_accept_state_result_return_void.hpp"
#include "rom/src/predicates/apply_jacobian_methods/rom_has_const_apply_jacobian_method_accept_state_operand_result_return_void.hpp"
#include "rom/src/predicates/preconditioning_methods/rom_has_const_apply_preconditioner_method_accept_state_result_return_void.hpp"

#include "rom/src/will_be_concepts/rom_masker.hpp"
#include "rom/src/will_be_concepts/rom_preconditioner.hpp"

#include "rom/src/will_be_concepts/system/rom_steady_system.hpp"
#include "rom/src/will_be_concepts/system/rom_discrete_time_system_with_user_provided_apply_jacobian.hpp"
#include "rom/src/will_be_concepts/system/rom_continuous_time_system_without_user_provided_apply_jacobian.hpp"
#include "rom/src/will_be_concepts/system/rom_continuous_time_system_with_user_provided_apply_jacobian.hpp"
#include "rom/src/will_be_concepts/system/rom_continuous_time_system.hpp"

//----------------
// galerkin
//----------------
#include "rom/src/galerkin/rom_compose_and_create_galerkin.hpp"

//-----------------
// LSPG
//-----------------
#include "rom/src/lspg/rom_default_lspg.hpp"
#include "rom/src/lspg/rom_preconditioned_default_lspg.hpp"
#include "rom/src/lspg/rom_masked_lspg.hpp"

//----------------
// wls
//----------------
#include "rom/src/wls/rom_wls_types.hpp"
#include "rom/src/wls/rom_wls_jacobian_updating_tag.hpp"
#include "rom/src/wls/rom_wls_jacobians_container.hpp"
#include "rom/src/wls/rom_wls_preconditioners.hpp"

#include "rom/src/wls/predicates/rom_wls_is_legitimate_preconditioner_type.hpp"
#include "rom/src/wls/predicates/rom_wls_is_legitimate_jacobian_updating_tag.hpp"
#include "rom/src/wls/time_schemes/rom_wls_implicit_euler.hpp"
#include "rom/src/wls/time_schemes/rom_wls_bdf2.hpp"
#include "rom/src/wls/time_schemes/rom_wls_select_timescheme_helper.hpp"

#include "rom/src/wls/rom_wls_hessian_gradient_system_api.hpp"
#include "rom/src/wls/rom_wls_hessian_and_gradient_sequential_policy.hpp"

#endif
