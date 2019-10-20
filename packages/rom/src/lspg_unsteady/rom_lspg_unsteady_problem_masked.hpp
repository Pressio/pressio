/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_problem_masked.hpp
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

#ifndef ROM_LSPG_UNSTEADY_PROBLEM_TYPE_GENERATOR_MASKED_HPP_
#define ROM_LSPG_UNSTEADY_PROBLEM_TYPE_GENERATOR_MASKED_HPP_

#include "./impl_velocity_api/rom_lspg_unsteady_type_generator_common_velocity_api.hpp"

namespace pressio{ namespace rom{

// template <
//   typename fom_type,
//   ode::ImplicitEnum odeName,
//   typename decoder_type,
//   typename lspg_state_type,
//   typename ud_ops = void,
//   mpl::enable_if_t<
//     ::pressio::rom::meta::model_meets_velocity_api_for_unsteady_lspg<fom_type>::value
//     > * = nullptr
//   >
// struct MaskedLSPGUnsteadyTypeGenerator
// {

//   // assert here that fom_type supports the mask interface

//   using common_types_helper =
//     impl::LSPGUnsteadyCommonTypesVelocityAPI<fom_type,
// 					     decoder_type,
// 					     lspg_state_type,
// 					     odeName,
// 					     ud_ops>;

//   using fom_t			= typename common_types_helper::fom_t;
//   using scalar_t		= typename common_types_helper::scalar_t;
//   using fom_native_state_t	= typename common_types_helper::fom_native_state_t;
//   using fom_state_t		= typename common_types_helper::fom_state_t;
//   using fom_velocity_t		= typename common_types_helper::fom_velocity_t;
//   using lspg_state_t		= typename common_types_helper::lspg_state_t;
//   using lspg_residual_t		= typename common_types_helper::lspg_residual_t;
//   using decoder_t		= typename common_types_helper::decoder_t;
//   using decoder_jac_t		= typename common_types_helper::decoder_jac_t;
//   using lspg_matrix_t		= typename common_types_helper::lspg_matrix_t;
//   using fom_state_reconstr_t	= typename common_types_helper::fom_state_reconstr_t;
//   using fom_states_data		= typename common_types_helper::fom_states_data;
//   using ud_ops_t		= typename common_types_helper::ud_ops_t;

//   // policy for evaluating the rhs of the fom object (<false> for unsteady overload)
//   using fom_eval_velocity_policy_t	= ::pressio::rom::policy::EvaluateFomVelocityDefault<false>;

//   // policy for left multiplying the fom jacobian with decoder_jac_t
//   // possibly involving other stuff like explained above (<false> for unsteady overload)
//   using fom_apply_jac_policy_t	= ::pressio::rom::policy::ApplyFomJacobianDefault<false>;

//   // policy defining how to compute the LSPG time-discrete residual
//   using lspg_residual_policy_t =
//     rom::decorator::Masked<
//     rom::LSPGUnsteadyResidualPolicyVelocityApi<
//       lspg_residual_t, fom_states_data, fom_eval_velocity_policy_t, ud_ops
//       >
//     >;

//   // policy defining how to compute the LSPG time-discrete jacobian
//   using lspg_jacobian_policy_t	=
//     rom::decorator::Masked<
//     rom::LSPGUnsteadyJacobianPolicyVelocityApi<
//       fom_states_data, lspg_matrix_t, fom_apply_jac_policy_t, decoder_t, ud_ops
//       >
//     >;

//   // auxiliary stepper
//   using aux_stepper_t = typename ::pressio::rom::impl::auxStepperHelper<
//     odeName, lspg_state_type,
//     lspg_residual_t, lspg_matrix_t,
//     fom_type, lspg_residual_policy_t,
//     lspg_jacobian_policy_t, scalar_t>::type;

//   // stepper object type
//   using lspg_stepper_t = ::pressio::ode::ImplicitStepper<
//     odeName, lspg_state_type,
//     lspg_residual_t, lspg_matrix_t,
//     fom_type, aux_stepper_t,
//     lspg_residual_policy_t, lspg_jacobian_policy_t, scalar_t>;

// };//end class


}}//end  namespace pressio::rom
#endif
