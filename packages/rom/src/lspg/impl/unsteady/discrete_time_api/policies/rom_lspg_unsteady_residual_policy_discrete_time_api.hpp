/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_residual_policy_discrete_time_api.hpp
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

#ifndef ROM_LSPG_IMPL_UNSTEADY_DISCRETE_TIME_API_POLICIES_ROM_LSPG_UNSTEADY_RESIDUAL_POLICY_DISCRETE_TIME_API_HPP_
#define ROM_LSPG_IMPL_UNSTEADY_DISCRETE_TIME_API_POLICIES_ROM_LSPG_UNSTEADY_RESIDUAL_POLICY_DISCRETE_TIME_API_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{ namespace unsteady{

template <
  typename residual_type,
  typename fom_states_manager_t,
  typename ud_ops_type = void
  >
class ResidualPolicyDiscreteTimeApi
{

public:
  using residual_t = residual_type;
  using ud_ops_t = ud_ops_type;

public:
  ResidualPolicyDiscreteTimeApi() = delete;
  ResidualPolicyDiscreteTimeApi(const ResidualPolicyDiscreteTimeApi &) = default;
  ResidualPolicyDiscreteTimeApi & operator=(const ResidualPolicyDiscreteTimeApi &) = delete;
  ResidualPolicyDiscreteTimeApi(ResidualPolicyDiscreteTimeApi &&) = default;
  ResidualPolicyDiscreteTimeApi & operator=(ResidualPolicyDiscreteTimeApi &&) = delete;
  ~ResidualPolicyDiscreteTimeApi() = default;

  // 1. void ops
  template <
    typename _ud_ops_t = ud_ops_t,
    ::pressio::mpl::enable_if_t< std::is_void<_ud_ops_t>::value, int > = 0
    >
  ResidualPolicyDiscreteTimeApi(fom_states_manager_t & fomStatesMngr)
    : fomStatesMngr_(fomStatesMngr){}

  // 2. nonvoid ops
  template <
    typename _ud_ops_t = ud_ops_t,
    ::pressio::mpl::enable_if_t< !std::is_void<_ud_ops_t>::value, int > = 0
    >
  ResidualPolicyDiscreteTimeApi(fom_states_manager_t & fomStatesMngr,
				const _ud_ops_t & udOps)
    : fomStatesMngr_(fomStatesMngr), udOps_(&udOps)
  {}

public:
  template <typename fom_system_t>
  residual_t create(const fom_system_t & system) const
  {
    residual_t R(system.createDiscreteTimeResidual());
    return R;
  }

  template <
    typename ode_tag,
    typename lspg_state_t,
    typename lspg_prev_states_t,
    typename fom_system_t,
    typename scalar_t
    >
  void compute(const lspg_state_t & romState,
	       const lspg_prev_states_t & romPrevStates,
	       const fom_system_t	& fomSystemObj,
	       const scalar_t		& time,
	       const scalar_t		& dt,
	       const ::pressio::ode::types::step_t & step,
	       residual_t & romR) const
  {
    this->compute_impl(romState, romPrevStates, fomSystemObj,
		       time, dt, step, romR);
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
    fomStatesMngr_.get().reconstructCurrentFomState(romState);

    /* the previous FOM states should only be recomputed when the time step changes
     * we do not need to reconstruct all the FOM states, we just need to reconstruct
     * the state at the previous step (i.e. t-dt) which is stored in romPrevStates[0]
     */
    if (storedStep_ != step){
      fomStatesMngr_.get() << romPrevStates.stateAt(ode::nMinusOne());
      storedStep_ = step;
    }
  }

  // we have here n = 1 prev rom states
  template <typename lspg_state_t, typename lspg_prev_states_t, typename fom_system_t, typename scalar_t>
  mpl::enable_if_t< lspg_prev_states_t::size()==1 >
  compute_impl(const lspg_state_t & romState,
	       const lspg_prev_states_t & romPrevStates,
	       const fom_system_t  & fomSystemObj,
	       const scalar_t & time,
	       const scalar_t & dt,
	       const ::pressio::ode::types::step_t & step,
	       residual_t & romR) const
  {
    doFomStatesReconstruction(romState, romPrevStates, step);
    const auto & yn   = fomStatesMngr_.get().currentFomStateCRef();
    const auto & ynm1 = fomStatesMngr_.get().fomStatePrevStepCRef();

    try{
      fomSystemObj.discreteTimeResidual(step, time, dt,
					*romR.data(),
					*yn.data(),
					*ynm1.data());
    }
    catch (::pressio::eh::discrete_time_residual_failure_unrecoverable const & e){
      throw ::pressio::eh::residual_evaluation_failure_unrecoverable();
    }
  }

  // we have here n = 2 prev rom states
  template <typename lspg_state_t, typename lspg_prev_states_t, typename fom_system_t, typename scalar_t>
  mpl::enable_if_t< lspg_prev_states_t::size()==2 >
  compute_impl(const lspg_state_t & romState,
	       const lspg_prev_states_t & romPrevStates,
	       const fom_system_t & fomSystemObj,
	       const scalar_t & time,
	       const scalar_t & dt,
	       const ::pressio::ode::types::step_t & step,
	       residual_t & romR) const
  {
    doFomStatesReconstruction(romState, romPrevStates, step);

    const auto & yn   = fomStatesMngr_.get().currentFomStateCRef();
    const auto & ynm1 = fomStatesMngr_.get().fomStatePrevStepCRef();
    const auto & ynm2 = fomStatesMngr_.get().fomStatePrevPrevStepCRef();

    try{
      fomSystemObj.discreteTimeResidual(step, time, dt,
					*romR.data(),
					*yn.data(),
					*ynm1.data(),
					*ynm2.data());

    }
    catch (::pressio::eh::discrete_time_residual_failure_unrecoverable const & e){
      throw ::pressio::eh::residual_evaluation_failure_unrecoverable();
    }
  }

protected:
  // storedStep is used to keep track of which step we are doing.
  // This is used to decide whether we need to update/recompute the previous
  // FOM states or not. Since it does not make sense to recompute previous
  // FOM states if we are not in a new time step.
  mutable ::pressio::ode::types::step_t storedStep_ = {};

  std::reference_wrapper<fom_states_manager_t> fomStatesMngr_;
  const ud_ops_t * udOps_ = {};
};

}}}}}//end namespace pressio::rom::lspg::unsteady::impl
#endif  // ROM_LSPG_IMPL_UNSTEADY_DISCRETE_TIME_API_POLICIES_ROM_LSPG_UNSTEADY_RESIDUAL_POLICY_DISCRETE_TIME_API_HPP_
