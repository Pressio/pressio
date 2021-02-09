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
    typename stencil_states_t,
    typename fom_system_t,
    typename scalar_t
    >
  void compute(const lspg_state_t & romState,
	       const stencil_states_t & stencilStates,
	       const fom_system_t & fomSystemObj,
	       const scalar_t & timeAtNextStep,
	       const scalar_t & dt,
	       const ::pressio::ode::types::step_t & currentStepNumber,
	       residual_type & romR) const
  {
    // since this is for hyp-red, I need to make sure the sTosInfo
    // is of the same extent as the romR
    assert(sTosInfo_.get().extent(0) == romR.extent(0));

    this->compute_impl<stepper_tag>(romState, romR, stencilStates,
				    fomSystemObj, timeAtNextStep, dt, currentStepNumber);
  }

  template <
    typename stepper_tag,
    typename lspg_state_t,
    typename stencil_states_t,
    typename fom_system_t,
    typename scalar_t,
    typename stencil_velocities_t
    >
  mpl::enable_if_t<
    std::is_same<stepper_tag, ::pressio::ode::implicitmethods::CrankNicolson>::value
    >
  compute(const lspg_state_t & romState,
	  const stencil_states_t & stencilStates,
	  const fom_system_t & fomSystemObj,
	  const scalar_t & timeAtNextStep,
	  const scalar_t & dt,
	  const ::pressio::ode::types::step_t & currentStepNumber,
	  stencil_velocities_t & stencilVelocities,
	  residual_type & romR) const
  {

    assert(sTosInfo_.get().extent(0) == romR.extent(0));

    this->compute_cn_impl<stepper_tag>
      (romState, romR, stencilStates, fomSystemObj,
       timeAtNextStep, dt, currentStepNumber, stencilVelocities);
  }

private:
  template <
    typename stepper_tag,
    typename lspg_state_t,
    typename stencil_states_t,
    typename fom_system_t,
    typename scalar_t
  >
  void compute_impl(const lspg_state_t & romState,
		    residual_type & romR,
		    const stencil_states_t & stencilStates,
		    const fom_system_t  & fomSystemObj,
		    const scalar_t & timeAtNextStep,
		    const scalar_t & dt,
		    const ::pressio::ode::types::step_t & currentStepNumber) const
  {
    /* the currrent FOM has to be recomputed every time regardless of
     * whether the currentStepNumber changes since we might be inside a non-linear solve
     * where the time step does not change but this residual method
     * is called multiple times.
     */
    fomStatesMngr_.get().reconstructAt(romState, ::pressio::ode::nPlusOne());

    /* previous FOM states should only be recomputed when the time step changes.
     * The method below does not recompute all previous states, but only
     * recomputes the n-th state and updates/shifts back all the other
     * FOM states stored. */
    if (storedStep_ != currentStepNumber){
      fomStatesMngr_.get().reconstructWithStencilUpdate(stencilStates(ode::n()));
      storedStep_ = currentStepNumber;
    }

    const auto & fomState = fomStatesMngr_(::pressio::ode::nPlusOne());
    fomSystemObj.velocity(*fomState.data(), timeAtNextStep, *romR.data());

    ::pressio::rom::lspg::impl::unsteady::time_discrete_residual
	<stepper_tag>(fomStatesMngr_.get(), romR, dt, sTosInfo_.get());
  }

  template <
    typename stepper_tag,
    typename lspg_state_t,
    typename stencil_states_t,
    typename fom_system_t,
    typename scalar_t,
    typename stencil_velocities_t
  >
  void compute_cn_impl(const lspg_state_t & romState,
		       residual_type & romR,
		       const stencil_states_t & stencilStates,
		       const fom_system_t & fomSystemObj,
		       const scalar_t & t_np1,
		       const scalar_t & dt,
		       const ::pressio::ode::types::step_t & currentStepNumber,
		       // for CN, stencilVelocities holds f_n+1 and f_n
		       stencil_velocities_t & stencilVelocities) const
  {
    PRESSIOLOG_DEBUG("residual policy with compute_cn_impl");

    fomStatesMngr_.get().reconstructAt(romState, ::pressio::ode::nPlusOne());

    if (storedStep_ != currentStepNumber){
      fomStatesMngr_.get().reconstructWithStencilUpdate(stencilStates(ode::n()));
      storedStep_ = currentStepNumber;

      // if the step changed, I need to compute f(y_n, t_n)
      const auto tn = t_np1-dt;
      auto & f_n = stencilVelocities(::pressio::ode::n());
      const auto & fomState_n = fomStatesMngr_(::pressio::ode::n());
      fomSystemObj.velocity(*fomState_n.data(), tn, *f_n.data());
    }

    // always compute f(y_n+1, t_n+1)
    auto & f_np1 = stencilVelocities(::pressio::ode::nPlusOne());
    const auto & fomState_np1 = fomStatesMngr_(::pressio::ode::nPlusOne());
    fomSystemObj.velocity(*fomState_np1.data(), t_np1, *f_np1.data());

    ::pressio::rom::lspg::impl::unsteady::time_discrete_residual
	<stepper_tag>(fomStatesMngr_.get(), stencilVelocities, romR, dt, sTosInfo_.get());
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
