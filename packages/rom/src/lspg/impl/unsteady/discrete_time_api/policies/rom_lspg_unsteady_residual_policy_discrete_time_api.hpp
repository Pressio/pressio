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

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{ namespace unsteady{ 

template <typename residual_type, typename fom_states_manager_t>
class ResidualPolicyDiscreteTimeApi
{

public:
  static constexpr bool isResidualPolicy_ = true;
  using residual_t = residual_type;

public:
  ResidualPolicyDiscreteTimeApi() = delete;
  ~ResidualPolicyDiscreteTimeApi() = default;

  ResidualPolicyDiscreteTimeApi(fom_states_manager_t & fomStatesMngr)
    : fomStatesMngr_(fomStatesMngr){}

public:
  template <typename system_t>
  residual_t create(const system_t & system) const
  {
    residual_t R(system.createDiscreteTimeResidual());
    return R;
  }

  template <
    typename ode_tag,
    typename lspg_state_t,
    typename lspg_prev_states_t,
    typename system_t,
    typename scalar_t
  >
  void compute(const lspg_state_t			& romState,
  		  const lspg_prev_states_t		& romPrevStates,
  		  const system_t				& system,
		  const scalar_t			& time,
		  const scalar_t			& dt,
		  const ::pressio::ode::types::step_t	& step,
		  residual_t				& romR,
		  ::pressio::Norm		normKind,
		  scalar_t				& normValue) const
  {
    this->compute_impl(romState, romPrevStates, system, time, dt, step, romR, normKind, normValue);
    // if (normKind == ::pressio::Norm::L2)
    //   normValue = ::pressio::ops::norm2(romR);
    // else if (normKind == ::pressio::Norm::L1)
    //   normValue = ::pressio::ops::norm1(romR);
    // else
    //   throw std::runtime_error("Invalid norm kind for lspg unsteady residual policy");
  }

private:
  template <typename lspg_state_t, typename lspg_prev_states_t>
  void doFomStatesReconstruction(const lspg_state_t & romState,
				 const lspg_prev_states_t & romPrevStates,
				 const ::pressio::ode::types::step_t & step) const
  {
    /* the currrent FOM has to be recomputed every time regardless
     * of whether the step changes since we might be inside a non-linear solve
     * where the time step does not change but this residual method
     * is called multiple times.
     */
    fomStatesMngr_.reconstructCurrentFomState(romState);

    /* the previous FOM states should only be recomputed when the time step changes
     * we do not need to reconstruct all the FOM states, we just need to reconstruct
     * the state at the previous step (i.e. t-dt) which is stored in romPrevStates[0]
     */
    if (storedStep_ != step){
      fomStatesMngr_ << romPrevStates.get(ode::nMinusOne());
      storedStep_ = step;
    }
  }

  // we have here n = 1 prev rom states
  template <typename lspg_state_t, typename lspg_prev_states_t, typename system_t, typename scalar_t>
  mpl::enable_if_t< lspg_prev_states_t::size()==1 >
  compute_impl(const lspg_state_t & romState,
	       const lspg_prev_states_t & romPrevStates,
	       const system_t  & system,
	       const scalar_t & time,
	       const scalar_t & dt,
	       const ::pressio::ode::types::step_t & step,
	       residual_t & romR, 
        ::pressio::Norm   normKind,
        scalar_t & normValue) const
  {
    doFomStatesReconstruction(romState, romPrevStates, step);
    const auto & yn   = fomStatesMngr_.getCRefToCurrentFomState();
    const auto & ynm1 = fomStatesMngr_.getCRefToFomStatePrevStep();
    ::pressio::rom::queryFomDiscreteTimeResidual(yn, ynm1, system, time, dt, step, romR, normKind, normValue);
  }

  // we have here n = 2 prev rom states
  template <typename lspg_state_t, typename lspg_prev_states_t, typename system_t, typename scalar_t>
  mpl::enable_if_t< lspg_prev_states_t::size()==2 >
  compute_impl(const lspg_state_t		        & romState,
		    const lspg_prev_states_t		& romPrevStates,
		    const system_t			        & system,
		    const scalar_t		        & time,
		    const scalar_t			& dt,
		    const ::pressio::ode::types::step_t & step,
		    residual_t			        & romR,
        ::pressio::Norm   normKind,
        scalar_t & normValue) const
  {
    doFomStatesReconstruction(romState, romPrevStates, step);

    const auto & yn   = fomStatesMngr_.getCRefToCurrentFomState();
    const auto & ynm1 = fomStatesMngr_.getCRefToFomStatePrevStep();
    const auto & ynm2 = fomStatesMngr_.getCRefToFomStatePrevPrevStep();
    ::pressio::rom::queryFomDiscreteTimeResidual(yn, ynm1, ynm2, system, time, dt, step, romR, normKind, normValue);
  }

protected:
  // storedStep is used to keep track of which step we are doing.
  // This is used to decide whether we need to update/recompute the previous
  // FOM states or not. Since it does not make sense to recompute previous
  // FOM states if we are not in a new time step.
  mutable ::pressio::ode::types::step_t storedStep_ = {};

  fom_states_manager_t & fomStatesMngr_;
};

}}}}}//end namespace pressio::rom::lspg::unsteady::impl
#endif
