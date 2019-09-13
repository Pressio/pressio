/*
//@HEADER
// ************************************************************************
//
// rom_lspg_residual_policy.hpp
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

#ifndef ROM_LSPG_RESIDUAL_POLICY_HPP_
#define ROM_LSPG_RESIDUAL_POLICY_HPP_

#include "../../rom_fwd.hpp"
#include "rom_lspg_time_discrete_residual.hpp"
#include "../../../../ode/src/implicit/policies/base/ode_implicit_residual_policy_base.hpp"
#include "../../rom_data_fom_rhs.hpp"
#include "../../rom_data_fom_states.hpp"

namespace pressio{ namespace rom{


template <
  typename fom_states_data,
  typename fom_rhs_data,
  typename fom_eval_rhs_policy,
  typename ud_ops
  >
class LSPGResidualPolicy
  : public ode::policy::ImplicitResidualPolicyBase<
      LSPGResidualPolicy<fom_states_data,
			 fom_rhs_data,
			 fom_eval_rhs_policy,
			 ud_ops>>,
    protected fom_states_data,
    protected fom_rhs_data,
    protected fom_eval_rhs_policy{

protected:
  using this_t = LSPGResidualPolicy<fom_states_data,
				    fom_rhs_data,
				    fom_eval_rhs_policy,
				    ud_ops>;
  friend ode::policy::ImplicitResidualPolicyBase<this_t>;

  using fom_states_data::yFom_;
  using fom_states_data::yFomOld_;
  using fom_states_data::maxNstates_;
  using fom_rhs_data::fomRhs_;

public:
  static constexpr bool isResidualPolicy_ = true;
  using typename fom_rhs_data::fom_rhs_t;

#ifdef HAVE_PYBIND11
  typename std::conditional<
    mpl::is_same<ud_ops, pybind11::object>::value,
    ud_ops,
    const ud_ops *
    >::type udOps_ = {};
#else 
    const ud_ops * udOps_ = {};
#endif


public:
  LSPGResidualPolicy() = delete;
  ~LSPGResidualPolicy() = default;

  // this cnstr only enabled when udOps is void
  template <
    typename _ud_ops = ud_ops,
    mpl::enable_if_t<
      std::is_void<_ud_ops>::value
      > * = nullptr
    >
  LSPGResidualPolicy(const fom_states_data & fomStates,
		     const fom_rhs_data & fomResids,
		     const fom_eval_rhs_policy & fomEvalRhsFunctor)
    : fom_states_data(fomStates),
      fom_rhs_data(fomResids),
      fom_eval_rhs_policy(fomEvalRhsFunctor){
    static_assert( std::is_void<_ud_ops>::value, "");
  }

  // this cnstr only enabled when udOps is non-void and not python
  template <
    typename _ud_ops = ud_ops,
    mpl::enable_if_t<
      !std::is_void<_ud_ops>::value 
#ifdef HAVE_PYBIND11
      and mpl::not_same<_ud_ops, pybind11::object>::value
#endif      
      > * = nullptr
    >
  LSPGResidualPolicy(const fom_states_data & fomStates,
  		     const fom_rhs_data & fomResids,
  		     const fom_eval_rhs_policy & fomEvalRhsFunctor,
  		     const _ud_ops & udOps)
    : fom_states_data(fomStates),
      fom_rhs_data(fomResids),
      fom_eval_rhs_policy(fomEvalRhsFunctor),
      udOps_{&udOps}{
    static_assert( !std::is_void<_ud_ops>::value, "");
  }

#ifdef HAVE_PYBIND11
  // this cnstr only enabled when udOps is non-void and python
  template <
    typename _ud_ops = ud_ops,
    mpl::enable_if_t<
      mpl::is_same<_ud_ops, pybind11::object>::value
      > * = nullptr
    >
  LSPGResidualPolicy(const fom_states_data & fomStates,
  		     const fom_rhs_data & fomResids,
  		     const fom_eval_rhs_policy & fomEvalRhsFunctor,
  		     const _ud_ops & udOps)
    : fom_states_data(fomStates),
      fom_rhs_data(fomResids),
      fom_eval_rhs_policy(fomEvalRhsFunctor),
      udOps_{udOps}
  {}
#endif
  

public:
  template <ode::ImplicitEnum odeMethod,
	    int n,
	    typename lspg_state_t,
	    typename lspg_residual_t,
	    typename fom_t,
	    typename scalar_t>
  void operator()(const lspg_state_t		   & romY,
		  lspg_residual_t		   & romR,
  		  const std::array<lspg_state_t,n> & oldYs,
  		  const fom_t			   & app,
		  scalar_t			   t,
		  scalar_t			   dt) const
  {
    this->compute_impl<odeMethod, n>(romY, romR, oldYs, app, t, dt);
  }

  template <ode::ImplicitEnum odeMethod,
	    int n,
	    typename lspg_state_t,
	    typename fom_t,
	    typename scalar_t>
  fom_rhs_t operator()(const lspg_state_t		   & romY,
			 const std::array<lspg_state_t,n>  & oldYs,
			 const fom_t			   & app,
			 scalar_t			   t,
			 scalar_t			   dt) const
  {
    this->compute_impl<odeMethod, n>(romY, fomRhs_, oldYs, app, t, dt);
    return fomRhs_;
  }

private:

  template <
    ode::ImplicitEnum odeMethod,
    int n,
    typename state_t,
    typename residual_t,
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
    rom::impl::time_discrete_residual
    <odeMethod, maxNstates_>(yFom, yFomOld, romR, dt);
  }

  template <
    ode::ImplicitEnum odeMethod,
    int n,
    typename state_t,
    typename residual_t,
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
    rom::impl::time_discrete_residual
      <odeMethod, maxNstates_>(yFom, yFomOld, romR, dt, udOps_);
  }

  template <ode::ImplicitEnum odeMethod,
	    int n,
	    typename lspg_state_t,
	    typename lspg_residual_t,
	    typename fom_t,
	    typename scalar_t>
  void compute_impl(const lspg_state_t		   & romY,
		    lspg_residual_t		   & romR,
		    const std::array<lspg_state_t,n> & oldYs,
		    const fom_t			   & app,
		    scalar_t			   t,
		    scalar_t			   dt) const
  {
#ifdef HAVE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("lspg residual");
#endif

    fom_states_data::template reconstructCurrentFomState(romY);
    fom_states_data::template reconstructFomOldStates<n>(oldYs);

#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("fom eval rhs");
#endif
    fom_eval_rhs_policy::evaluate(app, yFom_, romR, t);

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("fom eval rhs");
#endif

#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("time discrete residual");
#endif
    this->time_discrete_dispatcher<odeMethod,maxNstates_>(yFom_, yFomOld_, romR, dt);
#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("time discrete residual");
#endif

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("lspg residual");
#endif
  }


};//end class

}}//end namespace pressio::rom
#endif
