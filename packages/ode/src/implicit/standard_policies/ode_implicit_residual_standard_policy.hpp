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

template<
  typename state_type,
  typename system_type,
  typename residual_type
  >
class ResidualStandardPolicy<
  state_type, system_type, residual_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::ode::meta::legitimate_implicit_state_type<state_type>::value and
    ::pressio::ode::meta::legitimate_implicit_residual_type<residual_type>::value and
    containers::meta::is_wrapper<state_type>::value and
    containers::meta::is_wrapper<residual_type>::value
    >
  >
{

public:
  ResidualStandardPolicy() = default;
  ~ResidualStandardPolicy() = default;

public:
  residual_type operator()(const state_type & odeCurrentState, const system_type & model) const
  {
    residual_type R(model.velocity(*odeCurrentState.data(), 0.));
    return R;
  }

  template <typename tag_name, typename prev_states_type, typename scalar_type>
  void operator()(const state_type & odeCurrentState,
	     const prev_states_type & odePrevStates,
	     const system_type & model,
	     const scalar_type & t,
	     const scalar_type & dt,
	     const types::step_t & step,
	     residual_type & R,
	     ::pressio::solvers::Norm normKind,
	     scalar_type & normValue) const
  {
    model.velocity(*odeCurrentState.data(), t, *R.data());
    ::pressio::ode::impl::time_discrete_residual<tag_name>(odeCurrentState, R, odePrevStates, dt);
  }

};//end class

}}}}//end namespace pressio::ode::implicitmethods::policy
#endif
