/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_residual_policy_velocity_api.hpp
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

#ifndef ROM_LSPG_UNSTEADY_RESIDUAL_POLICY_VELOCITY_api_HPP_
#define ROM_LSPG_UNSTEADY_RESIDUAL_POLICY_VELOCITY_api_HPP_

#include "../../rom_fwd.hpp"
#include "../../rom_static_container_fom_states.hpp"
#include "../../../../ode/src/implicit/policies/base/ode_implicit_residual_policy_base.hpp"
#include "rom_lspg_time_discrete_residual.hpp"

namespace pressio{ namespace rom{ namespace lspg{ namespace unsteady{ namespace impl{

template <
  typename residual_type,
  typename fom_states_cont_type,
  typename fom_velocity_eval_policy,
  typename ud_ops
  >
class ResidualPolicyVelocityApi
  : public ::pressio::ode::policy::ImplicitResidualPolicyBase<
      ResidualPolicyVelocityApi<residual_type,
			 fom_states_cont_type,
			 fom_velocity_eval_policy,
			 ud_ops>>,
    protected fom_velocity_eval_policy
{

public:
  using this_t = ResidualPolicyVelocityApi<residual_type,
				    fom_states_cont_type,
				    fom_velocity_eval_policy,
				    ud_ops>;
  friend ::pressio::ode::policy::ImplicitResidualPolicyBase<this_t>;

  static constexpr bool isResidualPolicy_ = true;
  using residual_t = residual_type;

public:
  ResidualPolicyVelocityApi() = delete;
  ~ResidualPolicyVelocityApi() = default;

  /* for constructing this we need to deal with a few cases
   * 1. regular c++ with void ops
   * 2. regular c++ with non-void ops
   * 3. python bindings with void ops
   * 4. python bindings with non-void ops
   */

  // 1. enable for regular c++ and void ops
  template <
    typename _residual_type = residual_type,
    typename _ud_ops = ud_ops,
    ::pressio::mpl::enable_if_t<
      std::is_void<_ud_ops>::value
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
      and !::pressio::containers::meta::is_array_pybind11<_residual_type>::value
#endif
      > * = nullptr
    >
  ResidualPolicyVelocityApi(const _residual_type & RIn,
					fom_states_cont_type & fomStatesIn,
					const fom_velocity_eval_policy & fomEvalVelocityFunctor)
    : fom_velocity_eval_policy(fomEvalVelocityFunctor),
      R_{RIn},
      fomStates_(fomStatesIn)
  {}


  // 2. enable for regular c++ and non-void ops
  template <
    typename _residual_type = residual_type,
    typename _ud_ops = ud_ops,
    ::pressio::mpl::enable_if_t<
      !std::is_void<_ud_ops>::value
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
      and !::pressio::containers::meta::is_array_pybind11<_residual_type>::value
      and ::pressio::mpl::not_same<_ud_ops, pybind11::object>::value
#endif
      > * = nullptr
    >
  ResidualPolicyVelocityApi(const _residual_type & RIn,
					fom_states_cont_type & fomStatesIn,
					const fom_velocity_eval_policy & fomEvalVelocityFunctor,
					const _ud_ops & udOps)
    : R_{RIn},
      fomStates_(fomStatesIn),
      fom_velocity_eval_policy(fomEvalVelocityFunctor),
      udOps_{&udOps}
  {}


#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  // 3. python bindings with void ops (which means we do the ops here)
  template <
    typename _residual_type = residual_type,
    typename _ud_ops = ud_ops,
    ::pressio::mpl::enable_if_t<
      ::pressio::containers::meta::is_array_pybind11<_residual_type>::value and
      std::is_void<_ud_ops>::value
      > * = nullptr
    >
  ResidualPolicyVelocityApi(const _residual_type & RIn,
					fom_states_cont_type & fomStatesIn,
					const fom_velocity_eval_policy & fomEvalVelocityFunctor)
    : fom_velocity_eval_policy(fomEvalVelocityFunctor),
      R_{{_residual_type(const_cast<_residual_type &>(RIn).request())}},
      fomStates_(fomStatesIn)
  {}


  // 4. python bindings with non-void ops (means the user passes object with ops)
  template <
    typename _residual_type = residual_type,
    typename _ud_ops = ud_ops,
    ::pressio::mpl::enable_if_t<
      ::pressio::containers::meta::is_array_pybind11<_residual_type>::value and
      ::pressio::mpl::is_same<_ud_ops, pybind11::object>::value
      > * = nullptr
    >
  ResidualPolicyVelocityApi(const _residual_type & RIn,
					fom_states_cont_type & fomStatesIn,
					const fom_velocity_eval_policy & fomEvalVelocityFunctor,
					const _ud_ops & udOps)
    : fom_velocity_eval_policy(fomEvalVelocityFunctor),
      R_{{_residual_type(const_cast<_residual_type &>(RIn).request())}},
      fomStates_(fomStatesIn),
      udOps_{udOps}
  {}
#endif


public:
  template <
    ::pressio::ode::ImplicitEnum odeStepperName,
    std::size_t n,
    typename lspg_state_t,
    typename fom_t,
    typename scalar_t
  >
  void operator()(const lspg_state_t		   & romState,
		  residual_t			   & romR,
  		  const ::pressio::ode::StatesContainer<lspg_state_t,n> & romPrevStates,
  		  const fom_t			   & app,
		  const scalar_t		   & t,
		  const scalar_t		   & dt,
		  const ::pressio::ode::types::step_t & step) const
  {
    this->compute_impl<odeStepperName, n>(romState, romR, romPrevStates, app, t, dt, step);
  }

  template <
    ::pressio::ode::ImplicitEnum odeStepperName,
    std::size_t n,
    typename lspg_state_t,
    typename fom_t,
    typename scalar_t
    >
  residual_t operator()(const lspg_state_t		   & romState,
			const ::pressio::ode::StatesContainer<lspg_state_t,n>  & romPrevStates,
			const fom_t			   & app,
			const scalar_t			   & t,
			const scalar_t			   & dt,
			const ::pressio::ode::types::step_t & step) const
  {
    this->compute_impl<odeStepperName, n>(romState, R_, romPrevStates, app, t, dt, step);
    return R_;
  }



private:

  template <
    ::pressio::ode::ImplicitEnum odeStepperName,
    typename fom_state_cont_type,
    typename scalar_t,
    typename _ud_ops = ud_ops,
    ::pressio::mpl::enable_if_t<
	std::is_void<_ud_ops>::value
      > * = nullptr
  >
  void time_discrete_dispatcher(const fom_state_cont_type	& fomStates,
				residual_t			& romR,
				const scalar_t			& dt) const{
    using namespace ::pressio::rom::lspg::unsteady::impl;
    time_discrete_residual<odeStepperName>(fomStates, romR, dt);
  }

  template <
    ::pressio::ode::ImplicitEnum odeStepperName,
    typename fom_state_cont_type,
    typename scalar_t,
    typename _ud_ops = ud_ops,
    ::pressio::mpl::enable_if_t<
      !std::is_void<_ud_ops>::value
      > * = nullptr
  >
  void time_discrete_dispatcher(const fom_state_cont_type	& fomStates,
  				residual_t			& romR,
  				const scalar_t			& dt) const{
    using namespace ::pressio::rom::lspg::unsteady::impl;
    time_discrete_residual<odeStepperName>(fomStates, romR, dt, udOps_);
  }


private:
  template <
    ::pressio::ode::ImplicitEnum odeStepperName,
    std::size_t n,
    typename lspg_state_t,
    typename fom_t,
    typename scalar_t
  >
  void compute_impl(const lspg_state_t		     & romState,
		    residual_t			     & romR,
		    const ::pressio::ode::StatesContainer<lspg_state_t,n> & romPrevStates,
		    const fom_t			     & app,
		    const scalar_t		     & t,
		    const scalar_t		     & dt,
		    const ::pressio::ode::types::step_t & step) const
  {
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("lspg residual");
#endif

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

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("fom eval rhs");
#endif
    fom_velocity_eval_policy::evaluate(app, fomStates_.getCRefToCurrentFomState(), romR, t);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("fom eval rhs");
    timer->start("time discrete residual");
#endif

    this->time_discrete_dispatcher<odeStepperName>(fomStates_, romR, dt);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("time discrete residual");
    timer->stop("lspg residual");
#endif
  }


protected:
  // currentStep is used to keep track of which step we are doing.
  // This is used to decide whether we need to update/recompute the previous
  // FOM states or not. Since it does not make sense to recompute previous
  // FOM states if we are not in a new time step.
  mutable ::pressio::ode::types::step_t currentStep_ = {};

  mutable residual_t R_ = {};
  fom_states_cont_type & fomStates_;

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  // here we do this conditional type because it seems when ud_ops= pybind11::object
  // it only works if we copy the object. Need to figure out if we can leave ptr in all cases.
  typename std::conditional<
    ::pressio::mpl::is_same<ud_ops, pybind11::object>::value, ud_ops,
    const ud_ops *
    >::type udOps_ = {};
#else
    const ud_ops * udOps_ = {};
#endif

};//end class

}}}}}//end namespace pressio::rom::lspg::unstedy::impl
#endif
