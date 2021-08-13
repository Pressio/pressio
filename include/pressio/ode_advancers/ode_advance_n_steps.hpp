/*
//@HEADER
// ************************************************************************
//
// ode_advance_n_steps_explicit.hpp
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

#ifndef ODE_INTEGRATORS_ODE_ADVANCE_N_STEPS_EXPLICIT_HPP_
#define ODE_INTEGRATORS_ODE_ADVANCE_N_STEPS_EXPLICIT_HPP_

#include "./impl/ode_advance_noop_collector.hpp"
#include "./impl/ode_advance_n_steps.hpp"

namespace pressio{ namespace ode{

// explicit stepper
template<class StepperType, class StateType, class TimeType>
mpl::enable_if_t<
  ::pressio::ode::explicitly_steppable<StepperType, StateType, TimeType>::value
  >
advance_n_steps(StepperType & stepper,
		StateType & odeStateInOut,
		const TimeType start_time,
		const TimeType dt,
		const ::pressio::ode::step_count_type num_steps)
{

  using advancer_t  = impl::IntegratorNStepsWithConstDt;
  using collector_t = ::pressio::ode::impl::NoOpCollector<TimeType, StateType>;
  collector_t collector;
  advancer_t::execute(stepper, num_steps, start_time, dt, odeStateInOut, collector);
}

// explicit stepper, collector
template<class StepperType, class StateType, class TimeType, class collector_type>
mpl::enable_if_t<
  ::pressio::ode::explicitly_steppable<StepperType, StateType, TimeType>::value
  >
advance_n_steps(StepperType & stepper,
	      StateType & odeStateInOut,
	      const TimeType start_time,
	      const TimeType dt,
	      const ::pressio::ode::step_count_type num_steps,
	      collector_type & collector)
{

  static_assert
    (::pressio::ode::collector<collector_type, TimeType, StateType>::value,
     "You are trying to call advance_n_steps with an explicit stepper \
and a collector, but the collector type you are using is not admissible. \
It does not meet the API of a valid collector. \
See requirements in ode_is_legitimate_collector.hpp");

  using advancer_t  = impl::IntegratorNStepsWithConstDt;
  advancer_t::execute(stepper, num_steps, start_time, dt, odeStateInOut, collector);
}


/*
  implicit stepper
*/
template<
  class StepperType,
  class StateType,
  class TimeType,
  class SolverType,
  class ...Args
  >
mpl::enable_if_t<
  ::pressio::ode::implicitly_steppable<StepperType, StateType, TimeType, SolverType>::value and
  ::pressio::ode::legitimate_solver_for_implicit_stepper<SolverType, StepperType, StateType>::value
  >
advance_n_steps(StepperType & stepper,
	      StateType & odeStateInOut,
	      const TimeType startTime,
	      const TimeType dt,
	      const ::pressio::ode::step_count_type numSteps,
	      SolverType & solver,
	      Args && ...solver_args)
{

  using advancer_t  = impl::IntegratorNStepsWithConstDt;
  using collector_t = ::pressio::ode::impl::NoOpCollector<TimeType, StateType>;
  collector_t collector;
  advancer_t::execute(stepper, numSteps, startTime, dt,
		      odeStateInOut, collector,
		      solver, std::forward<Args>(solver_args)...);
}

/*
  implicit stepper and collector
*/
template<
  class StepperType,
  class StateType,
  class TimeType,
  class CollectorType,
  class SolverType,
  class ...Args
  >
::pressio::mpl::enable_if_t<
  ::pressio::ode::implicitly_steppable<StepperType, StateType, TimeType, SolverType>::value and
  ::pressio::ode::legitimate_solver_for_implicit_stepper<SolverType, StepperType, StateType>::value  and
  ::pressio::ode::collector<CollectorType, TimeType, StateType>::value
  >
advance_n_steps(StepperType & stepper,
	      StateType & odeStateInOut,
	      const TimeType startTime,
	      const TimeType dt,
	      const ::pressio::ode::step_count_type numSteps,
	      CollectorType & collector,
	      SolverType & solver,
	      Args && ...solver_args)
{

  using advancer_t  = impl::IntegratorNStepsWithConstDt;
  advancer_t::execute(stepper, numSteps, startTime, dt, odeStateInOut,
		      collector,
		      solver, std::forward<Args>(solver_args)...);
}

/*
  implicit stepper, solution guesser
*/
template<
  class StepperType,
  class StateType,
  class TimeType,
  class GuessCallbackType,
  class SolverType,
  class ...Args
  >
::pressio::mpl::enable_if_t<
  ::pressio::ode::implicitly_steppable_with_guesser<StepperType, StateType, TimeType, SolverType, GuessCallbackType>::value and
  ::pressio::ode::legitimate_solver_for_implicit_stepper<SolverType, StepperType, StateType>::value  and
  ::pressio::ode::is_legitimate_guesser<GuessCallbackType, ::pressio::ode::step_count_type, TimeType, StateType>::value
  >
advance_n_steps(StepperType & stepper,
		StateType & odeStateInOut,
		const TimeType startTime,
		const TimeType dt,
		const ::pressio::ode::step_count_type numSteps,
		GuessCallbackType && guessCb,
		SolverType & solver,
		Args && ...solver_args)

{

  using advancer_t  = impl::IntegratorNStepsWithConstDt;
  using collector_t = ::pressio::ode::impl::NoOpCollector<TimeType, StateType>;
  collector_t collector;
  advancer_t::execute(stepper, numSteps, startTime, dt, odeStateInOut,
		      collector, std::forward<GuessCallbackType>(guessCb),
		      solver, std::forward<Args>(solver_args)...);
}

/*
  implicit stepper, collector, guesser
*/
template<
  class StepperType,
  class StateType,
  class TimeType,
  class CollectorType,
  class SolverType,
  class GuessCallbackType,
  class ...Args
  >
::pressio::mpl::enable_if_t<
  ::pressio::ode::implicitly_steppable_with_guesser<StepperType, StateType, TimeType, SolverType, GuessCallbackType>::value and
  ::pressio::ode::collector<CollectorType, TimeType, StateType>::value and
  ::pressio::ode::legitimate_solver_for_implicit_stepper<SolverType, StepperType, StateType>::value and
  ::pressio::ode::is_legitimate_guesser<GuessCallbackType, ::pressio::ode::step_count_type, TimeType, StateType>::value
  >
advance_n_steps(StepperType & stepper,
	      StateType & odeStateInOut,
	      const TimeType startTime,
	      const TimeType dt,
	      const ::pressio::ode::step_count_type numSteps,
	      CollectorType & collector,
	      GuessCallbackType && guessCb,
	      SolverType & solver,
	      Args && ...solver_args)
{

  using advancer_t  = impl::IntegratorNStepsWithConstDt;
  advancer_t::execute(stepper, numSteps, startTime, dt, odeStateInOut,
		      collector, std::forward<GuessCallbackType>(guessCb),
		      solver, std::forward<Args>(solver_args)...);
}

// implicit stepper, dt setter
template<
  class StepperType,
  class StateType,
  class TimeType,
  class SolverType,
  class StepSizeSetterT,
  class ...Args
  >
::pressio::mpl::enable_if_t<
  ::pressio::ode::implicitly_steppable<StepperType, StateType, TimeType, SolverType>::value and
  ::pressio::ode::time_step_size_manager<StepSizeSetterT, ::pressio::ode::step_count_type, TimeType>::value and
  ::pressio::ode::legitimate_solver_for_implicit_stepper<SolverType, StepperType, StateType>::value
  >
advance_n_steps(StepperType & stepper,
	      StateType & odeStateInOut,
	      const TimeType start_time,
	      const ::pressio::ode::step_count_type num_steps,
	      StepSizeSetterT && dtManager,
	      SolverType & solver,
	      Args && ...solver_args)
{

  using advancer_t  = impl::IntegratorNStepsWithTimeStepSizeSetter;
  using collector_t = ::pressio::ode::impl::NoOpCollector<TimeType, StateType>;
  collector_t collector;
  advancer_t::execute(stepper, num_steps, start_time, odeStateInOut,
		      std::forward<StepSizeSetterT>(dtManager),
		      collector, solver,
		      std::forward<Args>(solver_args)...);
}

// implicit stepper, collector, with dt setter
template<
  class StepperType,
  class StateType,
  class TimeType,
  class StepSizeSetterType,
  class CollectorType,
  class SolverType,
  class ...Args
  >
::pressio::mpl::enable_if_t<
  ::pressio::ode::implicitly_steppable<StepperType, StateType, TimeType, SolverType>::value and
  ::pressio::ode::time_step_size_manager<StepSizeSetterType, ::pressio::ode::step_count_type, TimeType>::value and
  ::pressio::ode::collector<CollectorType, TimeType, StateType>::value and
  ::pressio::ode::legitimate_solver_for_implicit_stepper<SolverType, StepperType, StateType>::value
  >
advance_n_steps(StepperType & stepper,
		StateType & odeStateInOut,
		const TimeType start_time,
		const ::pressio::ode::step_count_type num_steps,
		StepSizeSetterType && dtManager,
		CollectorType & collector,
		SolverType & solver,
		Args && ...solver_args)
{

  using advancer_t  = impl::IntegratorNStepsWithTimeStepSizeSetter;
  advancer_t::execute(stepper, num_steps, start_time, odeStateInOut,
		      std::forward<StepSizeSetterType>(dtManager),
		      collector, solver,
		      std::forward<Args>(solver_args)...);
}


}}//end namespace pressio::ode
#endif  // ODE_INTEGRATORS_ODE_ADVANCE_N_STEPS_EXPLICIT_HPP_
