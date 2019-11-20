/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_residual_policy_residual_api.hpp
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

#ifndef ROM_LSPG_UNSTEADY_RESIDUAL_POLICY_RESIDUAL_API_HPP_
#define ROM_LSPG_UNSTEADY_RESIDUAL_POLICY_RESIDUAL_API_HPP_

#include "../../rom_fwd.hpp"
#include "../../../../ode/src/implicit/policies/base/ode_implicit_residual_policy_base.hpp"
#include "../../rom_static_container_fom_states.hpp"

namespace pressio{ namespace rom{ namespace lspg{ namespace unsteady{ namespace impl{

template <
  typename residual_type,
  typename fom_states_data_type,
  typename fom_querier_policy
  >
class ResidualPolicyResidualApi
  : public ode::policy::ImplicitResidualPolicyBase<
  ResidualPolicyResidualApi<residual_type, fom_states_data_type, fom_querier_policy>
  >,
  protected fom_querier_policy
{

public:
  using this_t = ResidualPolicyResidualApi<residual_type,
						       fom_states_data_type,
						       fom_querier_policy>;
  friend ode::policy::ImplicitResidualPolicyBase<this_t>;

  static constexpr bool isResidualPolicy_ = true;
  using residual_t = residual_type;

public:
  ResidualPolicyResidualApi() = delete;
  ~ResidualPolicyResidualApi() = default;

  ResidualPolicyResidualApi(fom_states_data_type & fomStatesIn,
					const fom_querier_policy & fomQuerierFunctor)
    : fom_querier_policy(fomQuerierFunctor),
      fomStates_(fomStatesIn){}

public:
  template <
    std::size_t n,
    typename lspg_state_t,
    typename fom_t,
    typename scalar_t
  >
  void operator()(const lspg_state_t			& romState,
  		  const ::pressio::ode::StatesContainer<lspg_state_t,n>	& romPrevStates,
  		  const fom_t				& app,
		  const scalar_t			& time,
		  const scalar_t			& dt,
		  const ::pressio::ode::types::step_t	& step,
		  residual_t				& romR) const
  {
    this->compute_impl(romState, romPrevStates, app, time, dt, step, romR);
  }

  template <typename lspg_state_t, typename fom_t>
  residual_t operator()(const lspg_state_t		    & romState,
			const fom_t			    & app) const
  {
    // this method only called once at the beginning
    fomStates_.template reconstructCurrentFomState(romState);
    residual_t R( fom_querier_policy::evaluate(fomStates_.getCRefToCurrentFomState(), app) );
    return R;
  }


private:
  template <std::size_t n, typename lspg_state_t>
  void doFomStatesReconstruction(const lspg_state_t & romState,
				 const ::pressio::ode::StatesContainer<lspg_state_t, n> & romPrevStates,
				 const ::pressio::ode::types::step_t & step) const
  {
    /* the currrent FOM has to be recomputed every time regardless
     * of whether the step changes since we might be inside a non-linear solve
     * where the time step does not change but this residual method
     * is called multiple times.
     */
    fomStates_.reconstructCurrentFomState(romState);

    /* the previous FOM states should only be recomputed when the time step changes
     * we do not need to reconstruct all the FOM states, we just need to reconstruct
     * the state at the previous step (i.e. t-dt) which is stored in romPrevStates[0]
     */
    if (currentStep_ != step){
      fomStates_ << romPrevStates[0];
      currentStep_ = step;
    }
  }

  // we have here n = 1 prev rom states
  template <typename lspg_state_t, typename fom_t, typename scalar_t>
  void compute_impl(const lspg_state_t		        & romState,
		    const ::pressio::ode::StatesContainer<lspg_state_t, 1> & romPrevStates,
		    const fom_t			        & app,
		    const scalar_t		        & time,
		    const scalar_t			& dt,
		    const ::pressio::ode::types::step_t & step,
		    residual_t			        & romR) const
  {
    doFomStatesReconstruction<1>(romState, romPrevStates, step);

    const auto & yn   = fomStates_.getCRefToCurrentFomState();
    const auto & ynm1 = fomStates_.getCRefToFomStatePrevStep();
    fom_querier_policy::evaluate(yn, ynm1, app, time, dt, step, romR);
  }

  // we have here n = 2 prev rom states
  template <typename lspg_state_t, typename fom_t, typename scalar_t>
  void compute_impl(const lspg_state_t		        & romState,
		    const ::pressio::ode::StatesContainer<lspg_state_t, 2> & romPrevStates,
		    const fom_t			        & app,
		    const scalar_t		        & time,
		    const scalar_t			& dt,
		    const ::pressio::ode::types::step_t & step,
		    residual_t			        & romR) const
  {
    doFomStatesReconstruction<2>(romState, romPrevStates, step);

    const auto & yn   = fomStates_.getCRefToCurrentFomState();
    const auto & ynm1 = fomStates_.getCRefToFomStatePrevStep();
    const auto & ynm2 = fomStates_.getCRefToFomStatePrevPrevStep();
    fom_querier_policy::evaluate(yn, ynm1, ynm2, app, time, dt, step, romR);
  }

protected:
  // currentStep is used to keep track of which step we are doing.
  // This is used to decide whether we need to update/recompute the previous
  // FOM states or not. Since it does not make sense to recompute previous
  // FOM states if we are not in a new time step.
  mutable ::pressio::ode::types::step_t currentStep_ = {};

  fom_states_data_type & fomStates_;
};

}}}}}//end namespace pressio::rom::lspg::unsteady::impl
#endif
