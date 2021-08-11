/*
//@HEADER
// ************************************************************************
//
// ode_advance_to_target_time_implicit_arbitrary_step_size.hpp
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

#ifndef ODE_INTEGRATORS_ODE_ADVANCE_TO_TARGET_TIME_HPP_
#define ODE_INTEGRATORS_ODE_ADVANCE_TO_TARGET_TIME_HPP_

#include "./impl/ode_advance_to_target_time.hpp"

namespace pressio{ namespace ode{

template<
  class StepperType,
  class StateType,
  class TimeType,
  class StepSizeSetterType,
  class SolverType,
  class ...Args
  >
::pressio::mpl::enable_if_t<
  ::pressio::ode::implicitly_steppable<StepperType, StateType, TimeType, SolverType>::value and
  ::pressio::ode::time_step_size_manager<StepSizeSetterType, ::pressio::ode::step_count_type, TimeType>::value and
  ::pressio::ode::legitimate_solver_for_implicit_stepper<SolverType, StepperType, StateType>::value
  >
advance_to_target_time(StepperType	& stepper,
		       StateType		& odeStateInOut,
		       const TimeType	start_time,
		       const TimeType	final_time,
		       StepSizeSetterType	&& dtManager,
		       SolverType		& solver,
		       Args&& ... solver_args)
{
  static_assert(::pressio::ode::implicit_state<StateType>::value,
		"You are trying to call advance_to_target_time with an implicit stepper \
but the state type you are using is not admissible for implicit time-stepping.");

  using collector_t = ::pressio::ode::impl::DummyCollector<TimeType, StateType>;
  collector_t collector;
  impl::integrate_to_target_time_with_time_step_size_manager<false>
    (stepper, start_time, final_time, odeStateInOut,
     collector, std::forward<StepSizeSetterType>(dtManager),
     solver, std::forward<Args>(solver_args)...);
}

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
advance_to_target_time(StepperType	& stepper,
		       StateType		& odeStateInOut,
		       const TimeType	start_time,
		       const TimeType	final_time,
		       StepSizeSetterType	&& dtManager,
		       CollectorType	& collector,
		       SolverType		& solver,
		       Args&& ...solver_args)
{

  static_assert
    (::pressio::ode::implicit_state<StateType>::value,
     "You are trying to call advance_to_target_time with an implicit stepper \
but the state type you are using is not admissible for implicit time-stepping.");

  impl::integrate_to_target_time_with_time_step_size_manager<false>
    (stepper, start_time, final_time, odeStateInOut,
     collector, std::forward<StepSizeSetterType>(dtManager),
     solver, std::forward<Args>(solver_args)...);
}

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
advance_to_target_time_with_time_step_recovery(StepperType	& stepper,
					       StateType	& odeStateInOut,
					       const TimeType	start_time,
					       const TimeType	final_time,
					       StepSizeSetterType	&& dtManager,
					       CollectorType	& collector,
					       SolverType	& solver,
					       Args&& ... solver_args)
{

  static_assert
    (::pressio::ode::implicit_state<StateType>::value,
     "You are trying to call advanceToTargetTime with an implicit stepper \
but the state type you are using is not admissible for implicit time-stepping.");

  impl::integrate_to_target_time_with_time_step_size_manager<true>
    (stepper, start_time, final_time, odeStateInOut,
     collector, std::forward<StepSizeSetterType>(dtManager),
     solver, std::forward<Args>(solver_args)...);
}

}}//end namespace pressio::ode
#endif  // ODE_INTEGRATORS_ODE_ADVANCE_TO_TARGET_TIME_IMPLICIT_ARBITRARY_STEP_SIZE_HPP_
