/*
//@HEADER
// ************************************************************************
//
// ode_time_discrete_residual.hpp
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

#ifndef ODE_TIME_DISCRETE_RESIDUAL_IMPL_HPP_
#define ODE_TIME_DISCRETE_RESIDUAL_IMPL_HPP_

#include "../ode_ConfigDefs.hpp"
#include "ode_implicit_constants.hpp"

namespace pressio{ namespace ode{ namespace impl{

// in all functions below:
// - on input R contains the application RHS, i.e. if the system is
// expressed as dudt = f(x,u,...), then on input R contains f(t_n, y_n, ...)
// - on output, R contains the time-discrete residual

template <
  ode::ImplicitEnum method,
  typename state_type,
  typename residual_type,
  typename pre_states_type,
  typename scalar_type
  >
mpl::enable_if_t<method == ode::ImplicitEnum::Euler>
time_discrete_residual(const state_type	& odeCurrentState,
		       residual_type	& R,
		       const pre_states_type & prevStates,
		       const scalar_type	& dt)
{
  using nm1 = ode::nMinusOne;

  constexpr auto cn   = ::pressio::ode::constants::bdf1<scalar_type>::c_n_;
  constexpr auto cnm1 = ::pressio::ode::constants::bdf1<scalar_type>::c_nm1_;
  const auto cf	  = ::pressio::ode::constants::bdf1<scalar_type>::c_f_ * dt;
  // R = y_n - y_n-1 - dt*f()
  ::pressio::containers::ops::do_update(R, cf, odeCurrentState, cn, prevStates.template get<nm1>(), cnm1);
}


template <
  ode::ImplicitEnum method,
  typename state_type,
  typename residual_type,
  typename pre_states_type,
  typename scalar_type
  >
mpl::enable_if_t< method == ode::ImplicitEnum::BDF2>
time_discrete_residual(const state_type	& odeCurrentState,
		       residual_type	& R,
		       const pre_states_type & prevStates,
		       const scalar_type	& dt)
{
  using nm1 = ode::nMinusOne;
  using nm2 = ode::nMinusTwo;

  constexpr auto cn   = ::pressio::ode::constants::bdf2<scalar_type>::c_n_;
  constexpr auto cnm1 = ::pressio::ode::constants::bdf2<scalar_type>::c_nm1_;
  constexpr auto cnm2 = ::pressio::ode::constants::bdf2<scalar_type>::c_nm2_;
  const auto cf	  = ::pressio::ode::constants::bdf2<scalar_type>::c_f_ * dt;

  // compute: R = y_n - 4/3 * y_n-1 + 1/3 * y_n-2 - 2/3 * dt * f(y_n, t_n)
  // R contains already f(y_n,t_n) so we can just update R by doing
  // R = -dt*2/3*R + y_n -4/3*y_n-1 + 1/3*y_n-2
  ::pressio::containers::ops::do_update(R, cf,
					odeCurrentState, cn,
					prevStates.template get<nm1>(), cnm1,
					prevStates.template get<nm2>(), cnm2);
}

}}}//end namespace pressio::ode::impl
#endif
