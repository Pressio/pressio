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

#ifndef ODE_INTEGRATORS_advance_N_STEPS_IMPLICIT_CONSTANT_STEP_SIZE_HPP_
#define ODE_INTEGRATORS_advance_N_STEPS_IMPLICIT_CONSTANT_STEP_SIZE_HPP_

#include "./impl/ode_call_stepper_policy.hpp"
#include "./impl/ode_n_steps_integrators.hpp"

namespace pressio{ namespace ode{

// basic version
template<
  typename stepper_type,
  typename state_type,
  typename time_type,
  typename solver_type
>
mpl::enable_if_t<
  ::pressio::ode::concepts::implicitly_steppable<stepper_type, state_type, time_type, solver_type>::value
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  and !::pressio::containers::predicates::is_array_pybind<state_type>::value
#endif
>
advanceNSteps(stepper_type & stepper,
		     state_type		 & odeStateInOut,
		     const time_type	 startTime,
		     const time_type	 dt,
		     const types::step_t numSteps,
		     solver_type	 & solver)
{

  static_assert(::pressio::ode::concepts::implicit_state<state_type>::value,
		"You are trying to call advanceNSteps with an implicit stepper \
but the state type you are using is not admissible for implicit time-stepping.");

  using do_step_policy_t = impl::ImplicitDoStepBasic<solver_type>;
  using advancer_t	 = impl::IntegratorNStepsWithConstDt<do_step_policy_t>;
  advancer_t::execute(numSteps, startTime, dt, odeStateInOut, stepper, solver);
}



#ifdef PRESSIO_ENABLE_TPL_PYBIND11
/* for pybind, we cannot use:
 * implicitmethods::StepperBase<stepper_type> & stepper
 * because the stepper is passed from Python so it does not know which overload to use
 * and we get a type error since it sees the stepper as an python object
 */

template<
  typename stepper_type,
  typename state_type,
  typename time_type,
  typename solver_type
>
::pressio::mpl::enable_if_t<
  ::pressio::ode::concepts::implicitly_steppable<stepper_type, state_type, time_type, solver_type>::value and
  ::pressio::containers::predicates::is_array_pybind<state_type>::value
>
advanceNSteps(stepper_type & stepper,
		     state_type		 & odeStateInOut,
		     const time_type	 startTime,
		     const time_type	 dt,
		     const types::step_t numSteps,
		     solver_type	 & solver)
{

  static_assert(::pressio::ode::concepts::implicit_state<state_type>::value,
		"You are trying to call advanceNSteps with an implicit stepper \
but the state type you are using is not admissible for implicit time-stepping.");

  // here we want to view the odeStateInOut since we want to modify its data,
  // which is numpy array owned by the user inside their Python code.
  // upon exit of this function, the original odeStateInOut is changed since odeStateView only views it.
  ::pressio::containers::Vector<state_type> odeStateView(odeStateInOut, ::pressio::view());

  using do_step_policy_t = impl::ImplicitDoStepBasic<solver_type>;
  using advancer_t	 = impl::IntegratorNStepsWithConstDt<do_step_policy_t>;
  advancer_t::execute(numSteps, startTime, dt, odeStateView, stepper, solver);
}
#endif


// with collector
template<
  typename stepper_type,
  typename state_type,
  typename time_type,
  typename collector_type,
  typename solver_type
>
::pressio::mpl::enable_if_t<
  ::pressio::ode::concepts::implicitly_steppable<stepper_type, state_type, time_type, solver_type>::value and
  ::pressio::ode::concepts::collector<collector_type, time_type, state_type>::value
>
advanceNSteps(stepper_type & stepper,
		     state_type		 & odeStateInOut,
		     const time_type	 startTime,
		     const time_type	 dt,
		     const types::step_t numSteps,
		     collector_type	 & collector,
		     solver_type	 & solver)
{

  static_assert(::pressio::ode::concepts::implicit_state<state_type>::value,
		"You are trying to call advanceNSteps with an implicit stepper \
but the state type you are using is not admissible for implicit time-stepping.");

  using do_step_policy_t = impl::ImplicitDoStepBasic<solver_type>;
  using advancer_t	 = impl::IntegratorNStepsWithCollectorAndConstDt<collector_type, do_step_policy_t>;
  advancer_t::execute(numSteps, startTime, dt, odeStateInOut, collector, stepper, solver);
}


// with guesser
template<
  typename stepper_type,
  typename state_type,
  typename time_type,
  typename solver_type,
  typename guess_callback_t
>
::pressio::mpl::enable_if_t<
  ::pressio::ode::concepts::implicitly_steppable_with_guesser<stepper_type, state_type, 
          time_type, solver_type, guess_callback_t>::value and
  ::pressio::ode::concepts::is_legitimate_guesser<
    guess_callback_t, types::step_t, time_type, state_type>::value
>
advanceNSteps(stepper_type & stepper,
		     state_type		 & odeStateInOut,
		     const time_type	 startTime,
		     const time_type	 dt,
		     const types::step_t numSteps,
		     solver_type	 & solver,
		     guess_callback_t && guessCb)
{

  static_assert(::pressio::ode::concepts::implicit_state<state_type>::value,
		"You are trying to call advanceNSteps with an implicit stepper \
but the state type you are using is not admissible for implicit time-stepping.");

  using do_step_policy_t = impl::ImplicitDoStepWithGuesser<solver_type, guess_callback_t>;
  using advancer_t	 = impl::IntegratorNStepsWithConstDt<do_step_policy_t>;
  advancer_t::execute(numSteps, startTime, dt, odeStateInOut, stepper, solver,
		      std::forward<guess_callback_t>(guessCb));
}


// with guesser and collector
template<
  typename stepper_type,
  typename state_type,
  typename time_type,
  typename collector_type,
  typename solver_type,
  typename guess_callback_t
>
::pressio::mpl::enable_if_t<
  ::pressio::ode::concepts::implicitly_steppable_with_guesser<stepper_type, state_type, 
          time_type, solver_type, guess_callback_t>::value and
  ::pressio::ode::concepts::collector<collector_type, time_type, state_type>::value and
  ::pressio::ode::concepts::is_legitimate_guesser<guess_callback_t, types::step_t, 
          time_type, state_type>::value
>
advanceNSteps(stepper_type & stepper,
		     state_type			& odeStateInOut,
		     const time_type	        startTime,
		     const time_type		dt,
		     const types::step_t	numSteps,
		     collector_type		& collector,
		     solver_type		& solver,
		     guess_callback_t		&& guessCb)
{

  static_assert(::pressio::ode::concepts::implicit_state<state_type>::value,
		"You are trying to call advanceNSteps with an implicit stepper \
but the state type you are using is not admissible for implicit time-stepping.");

  using do_step_policy_t = impl::ImplicitDoStepWithGuesser<solver_type, guess_callback_t>;
  using advancer_t	 = impl::IntegratorNStepsWithCollectorAndConstDt<collector_type, do_step_policy_t>;
  advancer_t::execute(numSteps, startTime, dt, odeStateInOut, collector, stepper,
		      solver, std::forward<guess_callback_t>(guessCb));
}

// with guesser and collector passed in different order
template<
  typename stepper_type,
  typename state_type,
  typename time_type,
  typename collector_type,
  typename solver_type,
  typename guess_callback_t
>
::pressio::mpl::enable_if_t<
  ::pressio::ode::concepts::implicitly_steppable_with_guesser<stepper_type, state_type, 
          time_type, solver_type, guess_callback_t>::value and
  ::pressio::ode::concepts::collector<collector_type, time_type, state_type>::value and
  ::pressio::ode::concepts::is_legitimate_guesser<guess_callback_t, types::step_t, 
          time_type, state_type>::value
>
advanceNSteps(stepper_type & stepper,
		     state_type			& odeStateInOut,
		     const time_type	        startTime,
		     const time_type		dt,
		     const types::step_t	numSteps,
		     solver_type		& solver,
		     collector_type		& collector,
		     guess_callback_t		&& guessCb)
{

  static_assert(::pressio::ode::concepts::implicit_state<state_type>::value,
		"You are trying to call advanceNSteps with an implicit stepper \
but the state type you are using is not admissible for implicit time-stepping.");

  advanceNSteps(stepper, odeStateInOut, startTime,
		  dt, numSteps, collector, solver,
		  std::forward<guess_callback_t>(guessCb));
}

}}//end namespace pressio::ode
#endif
