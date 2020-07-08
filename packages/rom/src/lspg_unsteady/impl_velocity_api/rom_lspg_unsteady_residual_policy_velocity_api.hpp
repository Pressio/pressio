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

#include "rom_lspg_time_discrete_residual.hpp"

namespace pressio{ namespace rom{ namespace lspg{ namespace unsteady{ namespace impl{

template <
  typename residual_type,
  typename fom_states_manager_t,
  typename ud_ops
  >
class ResidualPolicyVelocityApi
{

public:
  static constexpr bool isResidualPolicy_ = true;
  using residual_t = residual_type;

public:
  ResidualPolicyVelocityApi() = delete;
  ~ResidualPolicyVelocityApi() = default;

  /* for constructing this we need to deal with a few cases
   * 1. void ops
   * 2. non-void ops
   */

  // 1. void ops
  template <
    typename _residual_type = residual_type,
    typename _ud_ops = ud_ops,
    ::pressio::mpl::enable_if_t< std::is_void<_ud_ops>::value, int > = 0
    >
  ResidualPolicyVelocityApi(const _residual_type & RIn,
			    fom_states_manager_t & fomStatesMngr)
    : R_{RIn}, fomStatesMngr_(fomStatesMngr)
  {}

  // 2. non-void ops
  template <
    typename _residual_type = residual_type,
    typename _ud_ops = ud_ops,
    ::pressio::mpl::enable_if_t<
      !std::is_void<_ud_ops>::value, int > = 0
    >
  ResidualPolicyVelocityApi(const _residual_type & RIn,
			    fom_states_manager_t & fomStatesMngr,
			    const _ud_ops & udOps)
    : R_{RIn}, fomStatesMngr_(fomStatesMngr), udOps_{&udOps}
  {}

public:
  template <typename fom_t>
  residual_t create(const fom_t & app) const
  {
    return R_;
  }

  template <
    typename stepper_tag,
    typename lspg_state_t,
    typename prev_states_t,
    typename fom_t,
    typename scalar_t
  >
  void compute(const lspg_state_t & romState,
  		  const prev_states_t & romPrevStates,
  		  const fom_t & app,
		  const scalar_t & t,
		  const scalar_t & dt,
		  const ::pressio::ode::types::step_t & step,
		  residual_t & romR,
		  ::pressio::Norm normKind,
		  scalar_t & normValue) const
  {
    this->compute_impl<stepper_tag>(romState, romR, romPrevStates, app,
				    t, dt, step, normKind, normValue);
  }

private:
  template <
    typename stepper_tag,
    typename fom_state_cont_type,
    typename scalar_t,
    typename _ud_ops = ud_ops
  >
  ::pressio::mpl::enable_if_t< std::is_void<_ud_ops>::value >
  time_discrete_dispatcher(const fom_state_cont_type & fomStates,
			   residual_t & romR,
			   const scalar_t & dt,
			   ::pressio::Norm normKind,
			   scalar_t & normValue) const
  {
    using namespace ::pressio::rom::lspg::unsteady::impl;
    time_discrete_residual<stepper_tag>(fomStates, romR, dt);

    if (normKind == ::pressio::Norm::L2)
      normValue = ::pressio::ops::norm2(romR);
    else if (normKind == ::pressio::Norm::L1)
      normValue = ::pressio::ops::norm1(romR);
    else
      throw std::runtime_error("Invalid norm kind for lspg unsteady residual policy");
  }

  template <
    typename stepper_tag,
    typename fom_state_cont_type,
    typename scalar_t,
    typename _ud_ops = ud_ops
  >
  ::pressio::mpl::enable_if_t< !std::is_void<_ud_ops>::value >
  time_discrete_dispatcher(const fom_state_cont_type & fomStates,
			   residual_t & romR,
			   const scalar_t & dt,
			   ::pressio::Norm normKind,
			   scalar_t & normValue) const
  {
    using namespace ::pressio::rom::lspg::unsteady::impl;
    time_discrete_residual<stepper_tag>(fomStates, romR, dt, udOps_);

    if (normKind == ::pressio::Norm::L2)
      normValue = udOps_->norm2(*romR.data());
    else if (normKind == ::pressio::Norm::L1)
      normValue = udOps_->norm1(*romR.data());
    else
      throw std::runtime_error("Invalid norm kind for lspg unsteady residual policy");
  }

  template <
    typename stepper_tag,
    typename lspg_state_t,
    typename prev_states_t,
    typename fom_t,
    typename scalar_t
  >
  void compute_impl(const lspg_state_t		     & romState,
		    residual_t			     & romR,
		    const prev_states_t		     & romPrevStates,
		    const fom_t			     & app,
		    const scalar_t		     & t,
		    const scalar_t		     & dt,
		    const ::pressio::ode::types::step_t & step,
		    ::pressio::Norm normKind,
		    scalar_t & normValue) const
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
    fomStatesMngr_.reconstructCurrentFomState(romState);

    /* the previous FOM states should only be recomputed when the time step changes
     * we do not need to reconstruct all the FOM states, we just need to reconstruct
     * the state at the previous step (i.e. t-dt) which is stored in romPrevStates(0)
     */
    if (storedStep_ != step){
      fomStatesMngr_ << romPrevStates.get(ode::nMinusOne());
      storedStep_ = step;
    }

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->start("fom eval rhs");
#endif
    ::pressio::rom::queryFomVelocityUnsteady(app, fomStatesMngr_.getCRefToCurrentFomState(), romR, t);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("fom eval rhs");
    timer->start("time discrete residual");
#endif

    this->time_discrete_dispatcher<stepper_tag>(fomStatesMngr_, romR, dt, normKind, normValue);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("time discrete residual");
    timer->stop("lspg residual");
#endif
  }

protected:
  // storedStep is used to keep track of which step we are doing.
  // This is used to decide whether we need to update/recompute the previous
  // FOM states or not. Since it does not make sense to recompute previous
  // FOM states if we are not in a new time step.
  mutable ::pressio::ode::types::step_t storedStep_ = {};

  mutable residual_t R_ = {};
  fom_states_manager_t & fomStatesMngr_;

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
