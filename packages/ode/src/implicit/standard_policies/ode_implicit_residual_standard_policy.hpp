/*
//@HEADER
// ************************************************************************
//
// ode_implicit_residual_standard_policy.hpp
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

#ifndef ODE_POLICIES_STANDARD_IMPLICIT_RESIDUAL_STANDARD_POLICY_HPP_
#define ODE_POLICIES_STANDARD_IMPLICIT_RESIDUAL_STANDARD_POLICY_HPP_

namespace pressio{ namespace ode{ namespace implicitmethods{ namespace policy{

template<typename state_type, typename system_type, typename residual_type, typename=void>
class ResidualStandardPolicy;


// ---------------------------------------------------------------
// partially specialize for when continuous_time_system_implicit_stepping
// ---------------------------------------------------------------
template<
  typename state_type,
  typename system_type,
  typename residual_type
  >
class ResidualStandardPolicy<
  state_type, system_type, residual_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::ode::concepts::implicit_state<state_type>::value and
    ::pressio::ode::concepts::implicit_residual<residual_type>::value and
    containers::meta::is_wrapper<state_type>::value and
    containers::meta::is_wrapper<residual_type>::value and 
    ::pressio::ode::concepts::continuous_time_system_implicit_stepping<system_type>::value
    >
  >
{

public:
  ResidualStandardPolicy() = default;
  ~ResidualStandardPolicy() = default;

public:
  residual_type create(const system_type & system) const
  {
    residual_type R(system.createVelocity());
    return R;
  }

  template <typename ode_tag, typename prev_states_mgr_type, typename scalar_type>
  void compute(const state_type & odeCurrentState,
	     const prev_states_mgr_type & prevStatesMgr,
	     const system_type & system,
	     const scalar_type & t,
	     const scalar_type & dt,
	     const types::step_t & step,
	     residual_type & R,
	     ::pressio::Norm normKind,
	     scalar_type & normValue) const
  {
    system.velocity(*odeCurrentState.data(), t, *R.data());
    ::pressio::ode::impl::discrete_time_residual(odeCurrentState, R, prevStatesMgr, dt, ode_tag());
  }
};//end class


// ---------------------------------------------------------------
// partially specialize for when discrete_time_system_implicit_stepping
// ---------------------------------------------------------------
template<
  typename state_type,
  typename system_type,
  typename residual_type
  >
class ResidualStandardPolicy<
  state_type, system_type, residual_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::ode::concepts::implicit_state<state_type>::value and
    ::pressio::ode::concepts::implicit_residual<residual_type>::value and
    containers::meta::is_wrapper<state_type>::value and
    containers::meta::is_wrapper<residual_type>::value and
    ::pressio::ode::concepts::discrete_time_system_implicit_stepping<system_type>::value
    >
  >
{

public:
  ResidualStandardPolicy() = default;
  ~ResidualStandardPolicy() = default;

public:
  residual_type create(const system_type & system) const
  {
    residual_type R(system.createDiscreteTimeResidual());
    return R;
  }

  template <typename ode_tag, typename prev_states_mgr_type, typename scalar_type>
  mpl::enable_if_t< prev_states_mgr_type::size()==1>
  compute(const state_type & odeCurrentState,
      const prev_states_mgr_type & prevStatesMgr,
      const system_type & system,
      const scalar_type & t,
      const scalar_type & dt,
      const types::step_t & step,
      residual_type & R,
      ::pressio::Norm normKind,
      scalar_type & normValue) const
  {
    const auto & ynm1 = prevStatesMgr.get(ode::nMinusOne());
    system.template discreteTimeResidual(step, t, dt,
          *R.data(), normKind, normValue,
          *odeCurrentState.data(), *ynm1.data());
  }

  //-------------------------------
  // specialize for n == 2
  //-------------------------------
  template <typename ode_tag, typename prev_states_mgr_type, typename scalar_type>
  mpl::enable_if_t< prev_states_mgr_type::size()==2>
  compute(const state_type & odeCurrentState,
      const prev_states_mgr_type & prevStatesMgr,
      const system_type & system,
      const scalar_type & t,
      const scalar_type & dt,
      const types::step_t & step,
      residual_type & R,
      ::pressio::Norm normKind,
      scalar_type & normValue) const
  {
    const auto & ynm1 = prevStatesMgr.get(ode::nMinusOne());
    const auto & ynm2 = prevStatesMgr.get(ode::nMinusTwo());

    system.template discreteTimeResidual(step, t, dt,
          *R.data(), normKind, normValue,
          *odeCurrentState.data(), *ynm1.data(),*ynm2.data());
  }
};//end class

}}}}//end namespace pressio::ode::implicitmethods::policy
#endif
