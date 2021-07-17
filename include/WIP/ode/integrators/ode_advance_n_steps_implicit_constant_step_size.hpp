/*
//@HEADER
// ************************************************************************
//
// ode_advance_n_steps_implicit_constant_step_size.hpp
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

#ifndef ODE_INTEGRATORS_ODE_ADVANCE_N_STEPS_IMPLICIT_CONSTANT_STEP_SIZE_HPP_
#define ODE_INTEGRATORS_ODE_ADVANCE_N_STEPS_IMPLICIT_CONSTANT_STEP_SIZE_HPP_

#include "./impl/ode_n_steps_integrators.hpp"

namespace pressio{ namespace ode{

/*
  stepper, state, startTime, dt, # steps, solver, solver_args
*/
template<
  typename stepper_type,
  typename state_type,
  typename time_type,
  typename solver_type,
  typename ...Args
  >
mpl::enable_if_t<
  ::pressio::ode::constraints::implicitly_steppable<
    stepper_type, state_type, time_type, solver_type>::value and
  ::pressio::ode::constraints::legitimate_solver_for_implicit_stepper<
    solver_type, stepper_type, state_type>::value
  >
advanceNSteps(stepper_type & stepper,
	      state_type & odeStateInOut,
	      const time_type startTime,
	      const time_type dt,
	      const types::step_t numSteps,
	      solver_type & solver,
	      Args && ...solver_args)
{

  static_assert
    (::pressio::ode::constraints::implicit_state<state_type>::value,
     "You are trying to call advanceNSteps with an implicit stepper \
but the state type you are using is not admissible for implicit time-stepping.");

  using advancer_t  = impl::IntegratorNStepsWithConstDt;
  using collector_t = ::pressio::ode::impl::DummyCollector<time_type, state_type>;
  collector_t collector;
  advancer_t::execute(stepper, numSteps, startTime, dt,
		      odeStateInOut, collector,
		      solver, std::forward<Args>(solver_args)...);
}

/*
  stepper, state, startTime, dt, # steps, collector, solver, solver_args
*/
template<
  typename stepper_type,
  typename state_type,
  typename time_type,
  typename collector_type,
  typename solver_type,
  typename ...Args
  >
::pressio::mpl::enable_if_t<
  ::pressio::ode::constraints::implicitly_steppable<
    stepper_type, state_type, time_type, solver_type>::value and
  ::pressio::ode::constraints::legitimate_solver_for_implicit_stepper<
    solver_type, stepper_type, state_type>::value  and
  ::pressio::ode::constraints::collector<
    collector_type, time_type, state_type>::value
  >
advanceNSteps(stepper_type & stepper,
	      state_type & odeStateInOut,
	      const time_type startTime,
	      const time_type dt,
	      const types::step_t numSteps,
	      collector_type & collector,
	      solver_type & solver,
	      Args && ...solver_args)
{

  static_assert(::pressio::ode::constraints::implicit_state<state_type>::value,
		"You are trying to call advanceNSteps with an implicit stepper \
but the state type you are using is not admissible for implicit time-stepping.");

  using advancer_t  = impl::IntegratorNStepsWithConstDt;
  advancer_t::execute(stepper, numSteps, startTime, dt, odeStateInOut,
		      collector,
		      solver, std::forward<Args>(solver_args)...);
}

/*
  stepper, state, startTime, dt, # steps, guesser, solver, solver_args
*/
template<
  typename stepper_type,
  typename state_type,
  typename time_type,
  typename guess_callback_t,
  typename solver_type,
  typename ...Args
  >
::pressio::mpl::enable_if_t<
  ::pressio::ode::constraints::implicitly_steppable_with_guesser<
    stepper_type, state_type, time_type, solver_type, guess_callback_t>::value and
  ::pressio::ode::constraints::legitimate_solver_for_implicit_stepper<
    solver_type, stepper_type, state_type>::value  and
  ::pressio::ode::constraints::is_legitimate_guesser<
    guess_callback_t, types::step_t, time_type, state_type>::value
  >
advanceNSteps(stepper_type & stepper,
	      state_type & odeStateInOut,
	      const time_type startTime,
	      const time_type dt,
	      const types::step_t numSteps,
	      guess_callback_t && guessCb,
	      solver_type & solver,
	      Args && ...solver_args)

{

  static_assert(::pressio::ode::constraints::implicit_state<state_type>::value,
		"You are trying to call advanceNSteps with an implicit stepper \
but the state type you are using is not admissible for implicit time-stepping.");

  using advancer_t  = impl::IntegratorNStepsWithConstDt;
  using collector_t = ::pressio::ode::impl::DummyCollector<time_type, state_type>;
  collector_t collector;
  advancer_t::execute(stepper, numSteps, startTime, dt, odeStateInOut,
		      collector, std::forward<guess_callback_t>(guessCb),
		      solver, std::forward<Args>(solver_args)...);
}

/*
  stepper, state, startTime, dt, # steps, collector, guesser, solver, solver_args
*/
template<
  typename stepper_type,
  typename state_type,
  typename time_type,
  typename collector_type,
  typename solver_type,
  typename guess_callback_t,
  typename ...Args
  >
::pressio::mpl::enable_if_t<
  ::pressio::ode::constraints::implicitly_steppable_with_guesser<
    stepper_type, state_type, time_type, solver_type, guess_callback_t>::value and
  ::pressio::ode::constraints::collector<collector_type, time_type, state_type>::value and
  ::pressio::ode::constraints::legitimate_solver_for_implicit_stepper<
    solver_type, stepper_type, state_type>::value and
  ::pressio::ode::constraints::is_legitimate_guesser<
    guess_callback_t, types::step_t, time_type, state_type>::value
  >
advanceNSteps(stepper_type & stepper,
	      state_type & odeStateInOut,
	      const time_type startTime,
	      const time_type dt,
	      const types::step_t numSteps,
	      collector_type & collector,
	      guess_callback_t && guessCb,
	      solver_type & solver,
	      Args && ...solver_args)
{

  static_assert(::pressio::ode::constraints::implicit_state<state_type>::value,
		"You are trying to call advanceNSteps with an implicit stepper \
but the state type you are using is not admissible for implicit time-stepping.");

  using advancer_t  = impl::IntegratorNStepsWithConstDt;
  advancer_t::execute(stepper, numSteps, startTime, dt, odeStateInOut,
		      collector, std::forward<guess_callback_t>(guessCb),
		      solver, std::forward<Args>(solver_args)...);
}

}}//end namespace pressio::ode
#endif  // ODE_INTEGRATORS_ODE_ADVANCE_N_STEPS_IMPLICIT_CONSTANT_STEP_SIZE_HPP_
