/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_jacobian_policy_residual_api.hpp
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

#ifndef PRESSIO_ROM_EXP_GALERKIN_JACOBIAN_POLICY_RESIDUAL_API_HPP_
#define PRESSIO_ROM_EXP_GALERKIN_JACOBIAN_POLICY_RESIDUAL_API_HPP_

namespace pressio{ namespace rom{ namespace experimental{ namespace galerkin{ namespace impl{

template<
  typename rom_jacobian_type,
  typename fom_apply_jacobian_ret_type,
  typename decoder_type,
  typename fom_states_data_type,
  typename fom_querier_policy
  >
class JacobianPolicyResidualApi
{
  using scalar_t = typename ::pressio::containers::details::traits<rom_jacobian_type>::scalar_t;

public:
  using rom_jacobian_t = rom_jacobian_type;

public:
  JacobianPolicyResidualApi() = delete;
  ~JacobianPolicyResidualApi() = default;

  template< typename app_t>
  JacobianPolicyResidualApi(fom_states_data_type & fomStates,
			    const fom_querier_policy & fomQuerier,
			    const decoder_type & decoder,
			    const app_t & appObj)
    : fomQuerier_(fomQuerier),
      fomStates_(fomStates),
      phi_(decoder.getReferenceToJacobian()),
      fomApplyJac_(fomQuerier_.evaluate(fomStates_.getCRefToCurrentFomState(), appObj, phi_))
  {}

public:
  template <typename rom_state_t, typename rom_prev_states_t, typename fom_t>
  void operator()(const rom_state_t			& romState,
		  const rom_prev_states_t		& romPrevStates,
  		  const fom_t				& app,
		  const scalar_t			& time,
		  const scalar_t			& dt,
		  const ::pressio::ode::types::step_t	& step,
		  rom_jacobian_t			& romJac) const
  {
    this->compute_impl(romState, romPrevStates, app, time, dt, step, romJac);
  }

  // this is only called once
  template <typename rom_state_t, typename fom_t>
  rom_jacobian_t operator()(const rom_state_t & romState,
			    const fom_t	      & app) const
  {
    fomStates_.template reconstructCurrentFomState(romState);
    fom_apply_jacobian_ret_type fomApplyJac(fomQuerier_.evaluate(fomStates_.getCRefToCurrentFomState(), app, phi_));

    const auto nRows = phi_.extent(1);
    const auto nCols = phi_.extent(1);
    rom_jacobian_t romJac(nRows, nCols);
    constexpr auto zero = ::pressio::utils::constants::zero<scalar_t>();
    constexpr auto one  = ::pressio::utils::constants::one<scalar_t>();
    ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(),
                            one, phi_, fomApplyJac, zero, romJac);
    return romJac;
  }

private:
  // we have here n = 1 prev rom states
  template<
    typename rom_state_t, typename rom_prev_states_t, typename fom_t,
    typename scalar_t, typename rom_jac_t,
    mpl::enable_if_t< rom_prev_states_t::size()==1 > * = nullptr
  >
  void compute_impl(const rom_state_t			& romState,
		    const rom_prev_states_t		& romPrevStates,
  		    const fom_t			        & app,
  		    const scalar_t			& time,
  		    const scalar_t			& dt,
		    const ::pressio::ode::types::step_t	& step,
		    rom_jac_t				& romJac) const
  {
    // here we assume that the previous states have already been reconstructd
    // by the residual policy. So we do not recompute the FOM state.
    fomStates_.reconstructCurrentFomState(romState);

    const auto & yn   = fomStates_.getCRefToCurrentFomState();
    const auto & ynm1 = fomStates_.getCRefToFomStatePrevStep();
    fomQuerier_.evaluate(yn, ynm1, app, time, dt, step, phi_, fomApplyJac_);

    constexpr auto zero = ::pressio::utils::constants::zero<scalar_t>();
    constexpr auto one  = ::pressio::utils::constants::one<scalar_t>();
    ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(),
                            one, phi_, fomApplyJac_, zero, romJac);
  }

  // we have here n = 2 prev rom states
  template<
    typename rom_state_t, typename rom_prev_states_t, typename fom_t,
    typename scalar_t, typename rom_jac_t,
    mpl::enable_if_t< rom_prev_states_t::size()==2 > * = nullptr
    >
  void compute_impl(const rom_state_t			& romState,
		    const rom_prev_states_t		& romPrevStates,
  		    const fom_t			        & app,
  		    const scalar_t			& time,
  		    const scalar_t			& dt,
		    const ::pressio::ode::types::step_t	& step,
		    rom_jac_t				& romJac) const
  {
    // here we assume that the current state has already been reconstructd
    // by the residual policy. So we do not recompute the FOM state.
    // Maybe we should find a way to ensure this is the case.
    fomStates_.reconstructCurrentFomState(romState);

    const auto & yn   = fomStates_.getCRefToCurrentFomState();
    const auto & ynm1 = fomStates_.getCRefToFomStatePrevStep();
    const auto & ynm2 = fomStates_.getCRefToFomStatePrevStep();
    fomQuerier_.evaluate(yn, ynm1, ynm2, app, time, dt, step, phi_, fomApplyJac_);

    constexpr auto zero = ::pressio::utils::constants::zero<scalar_t>();
    constexpr auto one  = ::pressio::utils::constants::one<scalar_t>();
    ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(),
                            one, phi_, fomApplyJac_, zero, romJac);
  }

private:
  const fom_querier_policy & fomQuerier_;
  fom_states_data_type & fomStates_;
  const typename decoder_type::jacobian_type & phi_;
  mutable fom_apply_jacobian_ret_type fomApplyJac_;
};

}}}}}//end namespace
#endif
