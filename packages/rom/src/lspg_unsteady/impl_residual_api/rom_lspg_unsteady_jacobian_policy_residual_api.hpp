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

#ifndef ROM_LSPG_UNSTEADY_JACOBIAN_POLICY_RESIDUAL_api_HPP_
#define ROM_LSPG_UNSTEADY_JACOBIAN_POLICY_RESIDUAL_api_HPP_

#include "../../rom_fwd.hpp"
#include "../../../../ode/src/implicit/policies/base/ode_implicit_residual_policy_base.hpp"
#include "../../rom_static_container_fom_states.hpp"

namespace pressio{ namespace rom{ namespace impl{

template<
  typename fom_states_data_type,
  typename apply_jac_return_type,
  typename fom_querier_policy,
  typename decoder_type
  >
class LSPGUnsteadyJacobianPolicyResidualApi
  : public ::pressio::ode::policy::JacobianPolicyBase<
	LSPGUnsteadyJacobianPolicyResidualApi<fom_states_data_type,
			   apply_jac_return_type,
			   fom_querier_policy,
			   decoder_type>>,
    protected fom_querier_policy
{

public:
  using this_t = LSPGUnsteadyJacobianPolicyResidualApi<fom_states_data_type,
				    apply_jac_return_type,
				    fom_querier_policy,
				    decoder_type>;

  friend ::pressio::ode::policy::JacobianPolicyBase<this_t>;

  static constexpr bool isResidualPolicy_ = false;
  using apply_jac_return_t = apply_jac_return_type;

public:
  LSPGUnsteadyJacobianPolicyResidualApi() = delete;
  ~LSPGUnsteadyJacobianPolicyResidualApi() = default;

  LSPGUnsteadyJacobianPolicyResidualApi(fom_states_data_type & fomStates,
					const fom_querier_policy & fomQuerierFunctor,
					const decoder_type & decoder)
    : fom_querier_policy(fomQuerierFunctor),
      decoderObj_(decoder),
      fomStates_(fomStates){}

public:
  template <
    int n,
    typename lspg_state_t,
    typename lspg_jac_t,
    typename fom_t,
    typename scalar_t
  >
  void operator()(const lspg_state_t			& romState,
		  const ::pressio::ode::StatesContainer<lspg_state_t, n>  & romPrevStates,
  		  const fom_t				& app,
		  const scalar_t			& time,
		  const scalar_t			& dt,
		  const ::pressio::ode::types::step_t	& step,
		  lspg_jac_t				& romJac) const
  {
    this->compute_impl<n>(romState, romPrevStates, app, time, dt, step, romJac);
  }

  template <typename lspg_state_t, typename fom_t>
  apply_jac_return_t operator()(const lspg_state_t			& romState,
				const fom_t				& app) const
  {
    fomStates_.template reconstructCurrentFomState(romState);
    const auto & phi = decoderObj_.getReferenceToJacobian();
    apply_jac_return_t romJac(fom_querier_policy::evaluate(fomStates_.getCRefToCurrentFomState(), app, phi));
    return romJac;
  }

private:
  template<
    int n,
    typename lspg_state_t,
    typename fom_t,
    typename scalar_t,
    typename lspg_jac_t
  >
  void compute_impl(const lspg_state_t			& romState,
		    const ::pressio::ode::StatesContainer<lspg_state_t,n> & romPrevStates,
  		    const fom_t			        & app,
  		    const scalar_t			& time,
  		    const scalar_t			& dt,
		    const ::pressio::ode::types::step_t	& step,
		    lspg_jac_t				& romJac) const
  {
    // todo: this is not needed if jacobian is called after resiudal
    // because residual takes care of reconstructing the fom state
    //    timer->start("reconstruct fom state");
    fomStates_.template reconstructCurrentFomState(romState);

    const auto & phi = decoderObj_.getReferenceToJacobian();
    fom_querier_policy::template evaluate<fom_states_data_type::size()>(fomStates_, app, time, dt, step, phi, romJac);
  }

protected:
  const decoder_type & decoderObj_	= {};
  fom_states_data_type & fomStates_;
};

}}}//end namespace pressio::rom::impl
#endif
