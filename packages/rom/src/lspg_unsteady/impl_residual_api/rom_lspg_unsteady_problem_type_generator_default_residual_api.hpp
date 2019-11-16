/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_problem_type_generator_default_residual_api.hpp
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

#ifndef ROM_LSPG_UNSTEADY_PROBLEM_TYPE_GENERATOR_DEFAULT_RESIDUAL_API_HPP_
#define ROM_LSPG_UNSTEADY_PROBLEM_TYPE_GENERATOR_DEFAULT_RESIDUAL_API_HPP_

#include "rom_lspg_unsteady_residual_policy_residual_api.hpp"
#include "rom_lspg_unsteady_jacobian_policy_residual_api.hpp"
#include "../../fom_querying_policies/rom_query_fom_time_discrete_residual_policy.hpp"
#include "../../fom_querying_policies/rom_query_fom_apply_time_discrete_jacobian_policy.hpp"
#include "rom_lspg_unsteady_type_generator_common_residual_api.hpp"
#include "../../../../ode/src/ode_fwd.hpp"

namespace pressio{ namespace rom{ namespace lspg{ namespace unsteady{ namespace impl{

template <
  ::pressio::ode::ImplicitEnum odeName,
  typename fom_type,
  typename lspg_state_type,
  typename ... Args
  >
struct DefaultLSPGUnsteadyTypeGeneratorResidualApi
{
  static_assert( odeName == ::pressio::ode::ImplicitEnum::Arbitrary,
		 "\nTo use unsteady LSPG with the residual api, \
you need to pass ode::ImplicitEnum::Arbitrary");

  /* here, the fom_type must satisfy the residual api */
  static_assert( ::pressio::rom::meta::model_meets_residual_api_for_unsteady_lspg<fom_type>::value,
		 "\nUsing DefaultLSPGUnsteadyTypeGeneratorResidualApi \
requires a fom adapter class that meets the residual api. \
However, the fom/adapter type you passed does not meet this api. \
Verify the fom/adapter class you are using.");

  // pick the common types holder
  using common_types_t = LSPGUnsteadyCommonTypesResidualApi<odeName, fom_type, lspg_state_type, Args...>;

  using fom_t			= typename common_types_t::fom_t;
  using scalar_t		= typename common_types_t::scalar_t;
  using fom_native_state_t	= typename common_types_t::fom_native_state_t;
  using fom_native_residual_t	= typename common_types_t::fom_native_residual_t;

  using fom_state_t		= typename common_types_t::fom_state_t;
  using lspg_state_t		= typename common_types_t::lspg_state_t;
  using lspg_residual_t		= typename common_types_t::lspg_residual_t;

  using decoder_t		= typename common_types_t::decoder_t;
  using decoder_jac_t		= typename common_types_t::decoder_jac_t;
  using lspg_matrix_t		= typename common_types_t::lspg_matrix_t;
  using fom_state_reconstr_t	= typename common_types_t::fom_state_reconstr_t;
  using fom_states_data		= typename common_types_t::fom_states_data;
  using ud_ops_t		= typename common_types_t::ud_ops_t;

  // policy for evaluating the rhs of the fom object (<false> for unsteady overload)
  using fom_residual_querier_policy_t	= ::pressio::rom::policy::QueryFomTimeDiscreteResidual;

  // policy for querying the J*phi from FOM
  using fom_apply_jac_policy_t	= ::pressio::rom::policy::QueryFomApplyTimeDiscreteJacobian;

  // policy to compute the LSPG time-discrete residual
  using lspg_residual_policy_t	= ::pressio::rom::lspg::unsteady::impl::LSPGUnsteadyResidualPolicyResidualApi<
    lspg_residual_t, fom_states_data, fom_residual_querier_policy_t>;

  // policy to compute the LSPG time-discrete jacobian
  using lspg_jacobian_policy_t	= ::pressio::rom::lspg::unsteady::impl::LSPGUnsteadyJacobianPolicyResidualApi<
    fom_states_data, lspg_matrix_t, fom_apply_jac_policy_t, decoder_t>;

  using stepper_order_t  = typename common_types_t::order_setter;
  using tot_n_setter_t   = typename common_types_t::tot_n_setter;

  // declare type of stepper object
  using lspg_stepper_t		= ::pressio::ode::ImplicitStepper<
    odeName, lspg_state_t, lspg_residual_t, lspg_matrix_t, fom_type,
    lspg_residual_policy_t, lspg_jacobian_policy_t, stepper_order_t, tot_n_setter_t>;

};//end class

}}}}}//end  namespace pressio::rom::lspg::unsteady::impl
#endif
