/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_fom_residual_policy.hpp
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

#ifndef ROM_GALERKIN_IMPL_DISCRETE_TIME_API_POLICIES_ROM_GALERKIN_FOM_RESIDUAL_POLICY_HPP_
#define ROM_GALERKIN_IMPL_DISCRETE_TIME_API_POLICIES_ROM_GALERKIN_FOM_RESIDUAL_POLICY_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template <typename fom_states_manager_t, typename fom_residual_type>
class FomResidualPolicyDiscreteTimeApi
{
public:
  using data_type = fom_residual_type;

private:
  // storedStep is used to keep track of which step we are doing.
  // This is used to decide whether we need to update/recompute the previous
  // FOM states or not. Since it does not make sense to recompute previous
  // FOM states if we are not in a new time step.
  mutable ::pressio::ode::types::step_t storedStep_ = {};

  std::reference_wrapper<fom_states_manager_t> fomStatesMngr_;
  mutable fom_residual_type fomResidual_;

public:
  FomResidualPolicyDiscreteTimeApi() = delete;
  FomResidualPolicyDiscreteTimeApi(const FomResidualPolicyDiscreteTimeApi &) = default;
  FomResidualPolicyDiscreteTimeApi & operator=(const FomResidualPolicyDiscreteTimeApi &) = delete;
  FomResidualPolicyDiscreteTimeApi(FomResidualPolicyDiscreteTimeApi &&) = default;
  FomResidualPolicyDiscreteTimeApi & operator=(FomResidualPolicyDiscreteTimeApi &&) = delete;
  ~FomResidualPolicyDiscreteTimeApi() = default;

  template<typename fom_system_t>
  FomResidualPolicyDiscreteTimeApi(const fom_system_t & fomSystemObj,
				fom_states_manager_t & fomStatesMngr)
    : fomStatesMngr_(fomStatesMngr),
      fomResidual_(fomSystemObj.createDiscreteTimeResidual())
  {}

public:
  const fom_residual_type & get() const{
    return fomResidual_;
  }

  // we have here n = 1 stencil state
  template <
    typename galerkin_state_t,
    typename galerkin_stencil_states_t,
    typename fom_system_t,
    typename scalar_t
    >
  mpl::enable_if_t< galerkin_stencil_states_t::size()==1 >
  compute(const galerkin_state_t & galerkinState,
	  const fom_system_t & fomSystemObj,
	  const scalar_t & timeAtNextStep,
	  const scalar_t & dt,
	  const ::pressio::ode::types::step_t & currentStepNumber,
	  const galerkin_stencil_states_t & galerkinStencilStates) const
  {
    this->doFomStatesReconstruction(galerkinState, galerkinStencilStates, currentStepNumber);

    const auto & ynp1 = fomStatesMngr_(::pressio::ode::nPlusOne());
    const auto & yn   = fomStatesMngr_(::pressio::ode::n());
    fomSystemObj.discreteTimeResidual(currentStepNumber, timeAtNextStep, dt,
				      *fomResidual_.data(), *ynp1.data(), *yn.data());
  }

  // we have here n = 2 stencil state
  template <
    typename galerkin_state_t,
    typename galerkin_stencil_states_t,
    typename fom_system_t,
    typename scalar_t
    >
  mpl::enable_if_t< galerkin_stencil_states_t::size()==2 >
  compute(const galerkin_state_t & galerkinState,
	  const fom_system_t & fomSystemObj,
	  const scalar_t & timeAtNextStep,
	  const scalar_t & dt,
	  const ::pressio::ode::types::step_t & currentStepNumber,
	  const galerkin_stencil_states_t & galerkinStencilStates) const
  {
    this->doFomStatesReconstruction(galerkinState, galerkinStencilStates, currentStepNumber);

    const auto & ynp1 = fomStatesMngr_(::pressio::ode::nPlusOne());
    const auto & yn   = fomStatesMngr_(::pressio::ode::n());
    const auto & ynm1 = fomStatesMngr_(::pressio::ode::nMinusOne());

    fomSystemObj.discreteTimeResidual(currentStepNumber, timeAtNextStep,
				      dt, *fomResidual_.data(),
				      *ynp1.data(), *yn.data(), *ynm1.data());
  }

private:
  template <typename galerkin_state_t, typename galerkin_stencil_states_t>
  void doFomStatesReconstruction(const galerkin_state_t & galerkinState,
				 const galerkin_stencil_states_t & galerkinStencilStates,
				 const ::pressio::ode::types::step_t & currentStepNumber) const
  {
    /* the FOM state corresponding to the new predicted state has to be
     * recomputed every time regardless of the time step chaning or not,
     *  since we might be inside a non-linear solve
     * where the time step does not change but this residual method
     * is called multiple times. */
    fomStatesMngr_.get().reconstructAt(galerkinState, ::pressio::ode::nPlusOne());

    /* previous FOM states should only be recomputed when the time step changes.
     * The method below does not recompute all previous states, but only
     * recomputes the n-th state and updates/shifts back all the other
     * FOM states stored. */
    if (storedStep_ != currentStepNumber){
      fomStatesMngr_.get().reconstructWithStencilUpdate(galerkinStencilStates(ode::n()));
      storedStep_ = currentStepNumber;
    }
  }
};

}}}}//end namespace
#endif  // ROM_GALERKIN_IMPL_DISCRETE_TIME_API_POLICIES_ROM_GALERKIN_FOM_RESIDUAL_POLICY_HPP_
