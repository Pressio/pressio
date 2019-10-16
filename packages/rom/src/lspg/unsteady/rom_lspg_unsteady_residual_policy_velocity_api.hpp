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

#ifndef ROM_LSPG_UNSTEADY_RESIDUAL_POLICY_VELOCITY_API_HPP_
#define ROM_LSPG_UNSTEADY_RESIDUAL_POLICY_VELOCITY_API_HPP_

#include "../../rom_fwd.hpp"
#include "../../../../ode/src/implicit/policies/base/ode_implicit_residual_policy_base.hpp"
#include "../../rom_data_fom_states.hpp"
#include "rom_lspg_time_discrete_residual.hpp"

namespace pressio{ namespace rom{

template <
  typename residual_type,
  typename fom_states_data_type,
  typename fom_velocity_eval_policy,
  typename ud_ops
  >
class LSPGUnsteadyResidualPolicyVelocityApi
  : public ode::policy::ImplicitResidualPolicyBase<
      LSPGUnsteadyResidualPolicyVelocityApi<residual_type,
			 fom_states_data_type,
			 fom_velocity_eval_policy,
			 ud_ops>>,
    protected fom_velocity_eval_policy
{

public:
  using this_t = LSPGUnsteadyResidualPolicyVelocityApi<residual_type,
				    fom_states_data_type,
				    fom_velocity_eval_policy,
				    ud_ops>;
  friend ode::policy::ImplicitResidualPolicyBase<this_t>;

  static constexpr bool isResidualPolicy_ = true;
  using residual_t = residual_type;

public:
  LSPGUnsteadyResidualPolicyVelocityApi() = delete;
  ~LSPGUnsteadyResidualPolicyVelocityApi() = default;

  // cnstr enabled when udOps is void
  template <
    typename _ud_ops = ud_ops,
    mpl::enable_if_t<
      std::is_void<_ud_ops>::value
      > * = nullptr
    >
  LSPGUnsteadyResidualPolicyVelocityApi(const residual_t & RIn,
					fom_states_data_type & fomStatesIn,
					const fom_velocity_eval_policy & fomEvalVelocityFunctor)
    : R_{RIn},
      fomStates_(fomStatesIn),
      fom_velocity_eval_policy(fomEvalVelocityFunctor){
    static_assert( std::is_void<_ud_ops>::value, "");
  }

//   // cnstr enabled when udOps is non-void and not python
//   template <
//     typename _ud_ops = ud_ops,
//     mpl::enable_if_t<
//       !std::is_void<_ud_ops>::value
// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
//       and mpl::not_same<_ud_ops, pybind11::object>::value
// #endif
//       > * = nullptr
//     >
//   LSPGUnsteadyResidualPolicyVelocityApi(const fom_states_data & fomStates,
//   		     const fom_rhs_data & fomResids,
//   		     const fom_velocity_eval_policy & fomEvalRhsFunctor,
//   		     const _ud_ops & udOps)
//     : fom_states_data(fomStates),
//       fom_rhs_data(fomResids),
//       fom_velocity_eval_policy(fomEvalRhsFunctor),
//       udOps_{&udOps}{
//     static_assert( !std::is_void<_ud_ops>::value, "");
//   }

// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
//   // cnstr enabled when udOps is non-void and python
//   template <
//     typename _ud_ops = ud_ops,
//     mpl::enable_if_t<
//       mpl::is_same<_ud_ops, pybind11::object>::value
//       > * = nullptr
//     >
//   LSPGUnsteadyResidualPolicyVelocityApi(const fom_states_data & fomStates,
//   		     const fom_rhs_data & fomResids,
//   		     const fom_velocity_eval_policy & fomEvalRhsFunctor,
//   		     const _ud_ops & udOps)
//     : fom_states_data(fomStates),
//       fom_rhs_data(fomResids),
//       fom_velocity_eval_policy(fomEvalRhsFunctor),
//       udOps_{udOps}
//   {}
// #endif

public:
  template <
    ode::ImplicitEnum odeMethod,
    int n,
    typename lspg_state_t,
    typename fom_t,
    typename scalar_t
  >
  void operator()(const lspg_state_t		   & romY,
		  residual_t			   & romR,
  		  const std::array<lspg_state_t,n> & romOldYs,
  		  const fom_t			   & app,
		  scalar_t			   t,
		  scalar_t			   dt) const
  {
    this->compute_impl<odeMethod, n>(romY, romR, romOldYs, app, t, dt);
  }

  template <
    ode::ImplicitEnum odeMethod,
    int n,
    typename lspg_state_t,
    typename fom_t,
    typename scalar_t
    >
  residual_t operator()(const lspg_state_t		   & romY,
			const std::array<lspg_state_t,n>  & romOldYs,
			const fom_t			   & app,
			scalar_t			   t,
			scalar_t			   dt) const
  {
    this->compute_impl<odeMethod, n>(romY, R_, romOldYs, app, t, dt);
    return R_;
  }


private:
  template <
    ode::ImplicitEnum odeMethod,
    int n,
    typename state_t,
    typename scalar_t,
    typename _ud_ops = ud_ops,
    mpl::enable_if_t<
	std::is_void<_ud_ops>::value
      > * = nullptr
  >
  void time_discrete_dispatcher(const state_t			& yFom,
				const std::array<state_t,n>	& yFomOld,
				residual_t			& romR,
				scalar_t			dt) const{
    using namespace ::pressio::rom::impl;
    time_discrete_residual<odeMethod,
			   fom_states_data_type::N_>(yFom, yFomOld, romR, dt);
  }

  template <
    ode::ImplicitEnum odeMethod,
    int n,
    typename state_t,
    typename scalar_t,
    typename _ud_ops = ud_ops,
    mpl::enable_if_t<
      !std::is_void<_ud_ops>::value
      > * = nullptr
  >
  void time_discrete_dispatcher(const state_t			& yFom,
  				const std::array<state_t,n>	& yFomOld,
  				residual_t			& romR,
  				scalar_t			dt) const{
    using namespace ::pressio::rom::impl;
    time_discrete_residual<odeMethod,
			   fom_states_data_type::N_>(yFom, yFomOld, romR, dt, udOps_);
  }

  template <
    ode::ImplicitEnum odeMethod,
    int n,
    typename lspg_state_t,
    typename fom_t,
    typename scalar_t
  >
  void compute_impl(const lspg_state_t		     & romY,
		    residual_t			     & romR,
		    const std::array<lspg_state_t,n> & romOldYs,
		    const fom_t			     & app,
		    scalar_t			     t,
		    scalar_t			     dt) const
  {
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("lspg residual");
#endif

    fomStates_.template reconstructCurrentFomState(romY);
    fomStates_.template reconstructFomOldStates<n>(romOldYs);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("fom eval rhs");
#endif
    fom_velocity_eval_policy::evaluate(app, fomStates_.getCRefToFomState(), romR, t);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("fom eval rhs");
    timer->start("time discrete residual");
#endif

    this->time_discrete_dispatcher<odeMethod, fom_states_data_type::N_>
      (fomStates_.getCRefToFomState(), fomStates_.getCRefToFomOldStates(), romR, dt);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("time discrete residual");
    timer->stop("lspg residual");
#endif
  }


protected:
  mutable residual_t R_ = {};
  fom_states_data_type & fomStates_;

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  typename std::conditional<
    mpl::is_same<ud_ops, pybind11::object>::value,
    ud_ops, const ud_ops *
    >::type udOps_ = {};
#else
    const ud_ops * udOps_ = {};
#endif

};//end class

}}//end namespace pressio::rom
#endif
