/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_residual_policy_continuous_time_api.hpp
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

#ifndef ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_POLICIES_ROM_LSPG_UNSTEADY_RESIDUAL_POLICY_CONTINUOUS_TIME_API_HPP_
#define ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_POLICIES_ROM_LSPG_UNSTEADY_RESIDUAL_POLICY_CONTINUOUS_TIME_API_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{ namespace unsteady{

template <
  typename residual_type,
  typename fom_states_manager_t,
  typename ud_ops_type
  >
class ResidualPolicyContinuousTimeApi
{
public:
  using data_type = residual_type;

public:
  ResidualPolicyContinuousTimeApi() = delete;
  ResidualPolicyContinuousTimeApi(const ResidualPolicyContinuousTimeApi &) = default;
  ResidualPolicyContinuousTimeApi & operator=(const ResidualPolicyContinuousTimeApi &) = delete;
  ResidualPolicyContinuousTimeApi(ResidualPolicyContinuousTimeApi &&) = default;
  ResidualPolicyContinuousTimeApi & operator=(ResidualPolicyContinuousTimeApi &&) = delete;
  ~ResidualPolicyContinuousTimeApi() = default;

  /* for constructing this we need to deal with a few cases
   * 1. void ops
   * 2. non-void ops
   */

  // 1. void ops
  template <
    typename _ud_ops_t = ud_ops_type,
    ::pressio::mpl::enable_if_t< std::is_void<_ud_ops_t>::value, int > = 0
    >
  ResidualPolicyContinuousTimeApi(fom_states_manager_t & fomStatesMngr)
    : fomStatesMngr_(fomStatesMngr)
  {}

  // 2. non-void ops
  template <
    typename _ud_ops_t = ud_ops_type,
    ::pressio::mpl::enable_if_t<!std::is_void<_ud_ops_t>::value, int > = 0
    >
  ResidualPolicyContinuousTimeApi(fom_states_manager_t & fomStatesMngr,
				  const _ud_ops_t & udOps)
    : fomStatesMngr_(fomStatesMngr), udOps_{&udOps}
  {}

public:
  template <typename fom_system_t>
  residual_type create(const fom_system_t & fomObj) const
  {
    return residual_type(fomObj.createVelocity());
  }

  template <
    typename stepper_tag,
    typename lspg_state_t,
    typename stencil_states_t,
    typename fom_system_t,
    typename scalar_t
    >
  mpl::enable_if_t<
    mpl::not_same<stepper_tag, ::pressio::ode::implicitmethods::CrankNicolson>::value
    >
  compute(const lspg_state_t & romState,
	  const stencil_states_t & stencilStates,
	  const fom_system_t & fomSystemObj,
	  const scalar_t & timeAtNextStep,
	  const scalar_t & dt,
	  const ::pressio::ode::types::step_t & currentStepNumber,
	  residual_type & romR) const
  {
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
    this->compute_cn_impl<stepper_tag>
      (romState, romR, stencilStates, fomSystemObj,
       timeAtNextStep, dt, currentStepNumber, stencilVelocities);
  }

private:
  template <
  typename stepper_tag,
  typename fom_state_cont_type,
  typename scalar_t,
  typename _ud_ops_t = ud_ops_type
  >
  ::pressio::mpl::enable_if_t< std::is_void<_ud_ops_t>::value >
  time_discrete_dispatch(const fom_state_cont_type & fomStates,
			 residual_type & romR,
			 const scalar_t & dt) const
  {
    using namespace ::pressio::rom::lspg::impl::unsteady;
    time_discrete_residual<stepper_tag>(fomStates, romR, dt);
  }

  template <
    typename stepper_tag,
    typename fom_state_cont_type,
    typename scalar_t,
    typename _ud_ops_t = ud_ops_type
  >
  ::pressio::mpl::enable_if_t< !std::is_void<_ud_ops_t>::value >
  time_discrete_dispatch(const fom_state_cont_type & fomStates,
			 residual_type & romR,
			 const scalar_t & dt) const
  {
    using namespace ::pressio::rom::lspg::impl::unsteady;
    time_discrete_residual<stepper_tag>(fomStates, romR, dt, udOps_);
  }

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
		    const fom_system_t & fomSystemObj,
		    const scalar_t & timeAtNextStep,
		    const scalar_t & dt,
		    const ::pressio::ode::types::step_t & currentStepNumber) const
  {
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("lspg residual");
#endif

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

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("fom eval rhs");
#endif
    const auto & fomState = fomStatesMngr_(::pressio::ode::nPlusOne());
    fomSystemObj.velocity(*fomState.data(), timeAtNextStep, *romR.data());

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("fom eval rhs");
    timer->start("time discrete residual");
#endif

    this->time_discrete_dispatch<stepper_tag>(fomStatesMngr_.get(), romR, dt);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("time discrete residual");
    timer->stop("lspg residual");
#endif
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
	<stepper_tag>(fomStatesMngr_.get(), stencilVelocities, romR, dt);
  }

protected:
  // storedStep is used to keep track of which step we are doing.
  // This is used to decide whether we need to update/recompute the previous
  // FOM states or not. Since it does not make sense to recompute previous
  // FOM states if we are not in a new time step.
  mutable ::pressio::ode::types::step_t storedStep_ = {};

  std::reference_wrapper<fom_states_manager_t> fomStatesMngr_;

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  // here we do this conditional type because it seems when
  // ud_ops_t= pybind11::object it only works if we copy the object.
  // Need to figure out if we can leave ptr in all cases.
  typename std::conditional<
    ::pressio::mpl::is_same<ud_ops_type, pybind11::object>::value,
    ud_ops_type, const ud_ops_type *
    >::type udOps_ = {};
#else
  const ud_ops_type * udOps_ = {};
#endif

};//end class

}}}}}
#endif  // ROM_LSPG_IMPL_UNSTEADY_CONTINUOUS_TIME_API_POLICIES_ROM_LSPG_UNSTEADY_RESIDUAL_POLICY_CONTINUOUS_TIME_API_HPP_
