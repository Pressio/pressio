/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_hyper_reduced_residual_policy_continuous_time_api.hpp
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

#ifndef ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_POLICIES_ROM_LSPG_UNSTEADY_HYPER_REDUCED_RESIDUAL_POLICY_CONTINUOUS_TIME_API_HPP_
#define ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_POLICIES_ROM_LSPG_UNSTEADY_HYPER_REDUCED_RESIDUAL_POLICY_CONTINUOUS_TIME_API_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{ namespace unsteady{

template <
  typename residual_type,
  typename fom_states_manager_t,
  typename sample_to_stencil_t
  >
class HypRedResidualPolicyContinuousTimeApi
{

public:
  using data_type = residual_type;

public:
  HypRedResidualPolicyContinuousTimeApi() = delete;
  HypRedResidualPolicyContinuousTimeApi(const HypRedResidualPolicyContinuousTimeApi &) = default;
  HypRedResidualPolicyContinuousTimeApi & operator=(const HypRedResidualPolicyContinuousTimeApi &) = delete;
  HypRedResidualPolicyContinuousTimeApi(HypRedResidualPolicyContinuousTimeApi &&) = default;
  HypRedResidualPolicyContinuousTimeApi & operator=(HypRedResidualPolicyContinuousTimeApi &&) = delete;
  ~HypRedResidualPolicyContinuousTimeApi() = default;

  HypRedResidualPolicyContinuousTimeApi(fom_states_manager_t & fomStatesMngr,
					const sample_to_stencil_t & sTosInfo)
    : fomStatesMngr_(fomStatesMngr),
      sTosInfo_(sTosInfo)
  {}

public:
  template <typename fom_system_t>
  residual_type create(const fom_system_t & fomObj) const
  {
    return residual_type( fomObj.createVelocity() );
  }

  template <
    typename stepper_tag,
    typename lspg_state_t,
    typename prev_states_t,
    typename fom_system_t,
    typename scalar_t
    >
  void compute(const lspg_state_t & romState,
	       const prev_states_t & romPrevStates,
	       const fom_system_t & fomSystemObj,
	       const scalar_t & time,
	       const scalar_t & dt,
	       const ::pressio::ode::types::step_t & timeStep,
	       residual_type & romR) const
  {
    // since this is for hyp-red, I need to make sure the sTosInfo
    // is of the same extent as the romR
    assert(sTosInfo_.get().extent(0) == romR.extent(0));

    this->compute_impl<stepper_tag>(romState, romR, romPrevStates,
				    fomSystemObj, time, dt, timeStep);
  }

private:
  template <
    typename stepper_tag,
    typename lspg_state_t,
    typename prev_states_t,
    typename fom_system_t,
    typename scalar_t
  >
  void compute_impl(const lspg_state_t & romState,
		    residual_type & romR,
		    const prev_states_t & romPrevStates,
		    const fom_system_t  & fomSystemObj,
		    const scalar_t & time,
		    const scalar_t & dt,
		    const ::pressio::ode::types::step_t & timeStep) const
  {
    /* the currrent FOM has to be recomputed every time regardless of
     * whether the timeStep changes since we might be inside a non-linear solve
     * where the time step does not change but this residual method
     * is called multiple times.
     */
    fomStatesMngr_.get().reconstructCurrentFomState(romState);

    /* the previous FOM states should only be recomputed when time step changes.
     * no need to reconstruct all the FOM states, we just need to reconstruct
     * the state at the previous step (i.e. t-dt)
     */
    if (storedStep_ != timeStep){
      fomStatesMngr_.get() << romPrevStates.stateAt(ode::nMinusOne());
      storedStep_ = timeStep;
    }
    const auto & currentFomState = fomStatesMngr_.get().currentFomStateCRef();
    fomSystemObj.velocity(*currentFomState.data(), time, *romR.data());

    ::pressio::rom::lspg::impl::unsteady::time_discrete_residual
	<stepper_tag>(fomStatesMngr_.get(), romR, dt, sTosInfo_.get());
  }

protected:
  // storedStep is used to keep track of which step we are doing.
  // This is used to decide whether we need to update/recompute the previous
  // FOM states or not. Since it does not make sense to recompute previous
  // FOM states if we are not in a new time step.
  mutable ::pressio::ode::types::step_t storedStep_ = {};

  std::reference_wrapper<fom_states_manager_t> fomStatesMngr_;
  std::reference_wrapper<const sample_to_stencil_t> sTosInfo_;
};

}}}}}
#endif  // ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_POLICIES_ROM_LSPG_UNSTEADY_HYPER_REDUCED_RESIDUAL_POLICY_CONTINUOUS_TIME_API_HPP_
