/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_default_problem_traits_discrete_time_api.hpp
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

#ifndef PRESSIO_ROM_EXP_GALERKIN_PROBLEM_TYPE_GENERATOR_DEFAULT_RESIDUAL_API_IMPL_HPP_
#define PRESSIO_ROM_EXP_GALERKIN_PROBLEM_TYPE_GENERATOR_DEFAULT_RESIDUAL_API_IMPL_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template <
  typename stepper_tag,
  typename fom_type,
  typename rom_state_type,
  typename rom_jacobian_type,
  typename ... Args
  >
struct DefaultProblemTraitsDiscreteTimeApi
{
  // pick the common types holder
  using common_types = ::pressio::rom::galerkin::impl::CommonTraitsDiscreteTimeApi<
        fom_type, rom_state_type, rom_jacobian_type, Args...>;

  using fom_t			= typename common_types::fom_t;
  using scalar_t		= typename common_types::scalar_t;
  using fom_native_state_t	= typename common_types::fom_native_state_t;
  using fom_native_residual_t	= typename common_types::fom_native_residual_t;

  using fom_state_t		= typename common_types::fom_state_t;
  using fom_residual_t		= typename common_types::fom_residual_t;
  using fom_apply_jacobian_t	= typename common_types::fom_apply_jacobian_t;

  using rom_state_t		= typename common_types::rom_state_t;
  using rom_residual_t		= typename common_types::rom_residual_t;
  using rom_jacobian_t		= typename common_types::rom_jacobian_t;

  using decoder_t		= typename common_types::decoder_t;
  using decoder_jac_t		= typename common_types::decoder_jac_t;
  using fom_state_reconstr_t	= typename common_types::fom_state_reconstr_t;
  using fom_states_manager_t		= typename common_types::fom_states_manager_t;
  using ud_ops_t		= typename common_types::ud_ops_t;

  using residual_policy_t	= ::pressio::rom::galerkin::impl::ResidualPolicyDiscreteTimeApi<
    rom_residual_t, fom_residual_t, decoder_t, fom_states_manager_t>;

  using jacobian_policy_t	= ::pressio::rom::galerkin::impl::JacobianPolicyDiscreteTimeApi<
    rom_jacobian_t, fom_apply_jacobian_t, decoder_t, fom_states_manager_t>;

  using stepper_order_t  = typename common_types::order_setter;
  using tot_n_setter_t   = typename common_types::tot_n_setter;

  // stepper object
  using stepper_t = ::pressio::ode::ImplicitStepper<
    stepper_tag, rom_state_t, rom_residual_t, rom_jacobian_t, fom_type,
    stepper_order_t, tot_n_setter_t, residual_policy_t, jacobian_policy_t>;

};//end class

}}}}//end namespace
#endif
