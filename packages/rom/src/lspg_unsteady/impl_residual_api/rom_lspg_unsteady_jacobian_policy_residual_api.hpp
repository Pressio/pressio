/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_jacobian_policy_residual_api.hpp
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

#ifndef ROM_LSPG_UNSTEADY_JACOBIAN_POLICY_RESIDUAL_API_HPP_
#define ROM_LSPG_UNSTEADY_JACOBIAN_POLICY_RESIDUAL_API_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace unsteady{ namespace impl{

template<
  typename fom_states_data_type,
  typename apply_jac_return_type,
  typename fom_querier_policy,
  typename decoder_type
  >
class JacobianPolicyResidualApi
{

public:
  static constexpr bool isResidualPolicy_ = false;
  using apply_jac_return_t = apply_jac_return_type;

public:
  JacobianPolicyResidualApi() = delete;
  ~JacobianPolicyResidualApi() = default;

  JacobianPolicyResidualApi(fom_states_data_type & fomStates,
			    const fom_querier_policy & fomQuerier,
			    const decoder_type & decoder)
    : fomQuerier_(fomQuerier), decoderObj_(decoder), fomStates_(fomStates){}

public:
  template <
    typename lspg_state_t,
    typename lspg_prev_states_t,
    typename lspg_jac_t,
    typename fom_t,
    typename scalar_t
  >
  void operator()(const lspg_state_t			& romState,
		  const lspg_prev_states_t		& romPrevStates,
  		  const fom_t				& app,
		  const scalar_t			& time,
		  const scalar_t			& dt,
		  const ::pressio::ode::types::step_t	& step,
		  lspg_jac_t				& romJac) const
  {
    this->compute_impl(romState, romPrevStates, app, time, dt, step, romJac);
  }

  template <typename lspg_state_t, typename fom_t>
  apply_jac_return_t operator()(const lspg_state_t			& romState,
				const fom_t				& app) const
  {
    // this is only called once
    fomStates_.template reconstructCurrentFomState(romState);
    const auto & phi = decoderObj_.getReferenceToJacobian();
    apply_jac_return_t romJac(fomQuerier_.evaluate(fomStates_.getCRefToCurrentFomState(), app, phi));
    return romJac;
  }

private:

  // we have here n = 1 prev rom states
  template<
    typename lspg_state_t, typename lspg_prev_states_t, typename fom_t,
    typename scalar_t, typename lspg_jac_t,
    mpl::enable_if_t< lspg_prev_states_t::size()==1 > * = nullptr
  >
  void compute_impl(const lspg_state_t			& romState,
		    const lspg_prev_states_t		& romPrevStates,
  		    const fom_t			        & app,
  		    const scalar_t			& time,
  		    const scalar_t			& dt,
		    const ::pressio::ode::types::step_t	& step,
		    lspg_jac_t				& romJac) const
  {
    // here we assume that the current state has already been reconstructd
    // by the residual policy. So we do not recompute the FOM state.
    // Maybe we should find a way to ensure this is the case.
    fomStates_.template reconstructCurrentFomState(romState);

    const auto & phi = decoderObj_.getReferenceToJacobian();
    const auto & yn   = fomStates_.getCRefToCurrentFomState();
    const auto & ynm1 = fomStates_.getCRefToFomStatePrevStep();
    fomQuerier_.evaluate(yn, ynm1, app, time, dt, step, phi, romJac);
  }


  // we have here n = 2 prev rom states
  template<
    typename lspg_state_t, typename lspg_prev_states_t, typename fom_t,
    typename scalar_t, typename lspg_jac_t,
    mpl::enable_if_t< lspg_prev_states_t::size()==2 > * = nullptr
    >
  void compute_impl(const lspg_state_t			& romState,
		    const lspg_prev_states_t		& romPrevStates,
  		    const fom_t			        & app,
  		    const scalar_t			& time,
  		    const scalar_t			& dt,
		    const ::pressio::ode::types::step_t	& step,
		    lspg_jac_t				& romJac) const
  {
    // here we assume that the current state has already been reconstructd
    // by the residual policy. So we do not recompute the FOM state.
    // Maybe we should find a way to ensure this is the case.
    fomStates_.template reconstructCurrentFomState(romState);

    const auto & phi = decoderObj_.getReferenceToJacobian();
    const auto & yn   = fomStates_.getCRefToCurrentFomState();
    const auto & ynm1 = fomStates_.getCRefToFomStatePrevStep();
    const auto & ynm2 = fomStates_.getCRefToFomStatePrevStep();
    fomQuerier_.evaluate(yn, ynm1, ynm2, app, time, dt, step, phi, romJac);
  }

protected:
  const fom_querier_policy & fomQuerier_;
  const decoder_type & decoderObj_	= {};
  fom_states_data_type & fomStates_;
};

}}}}}//end namespace pressio::rom::lspg::unsteady::impl
#endif
