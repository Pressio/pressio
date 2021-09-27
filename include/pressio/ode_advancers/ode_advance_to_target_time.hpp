/*
//@HEADER
// ************************************************************************
//
// ode_advance_to_target_time.hpp
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

#ifndef ODE_ADVANCERS_ODE_ADVANCE_TO_TARGET_TIME_HPP_
#define ODE_ADVANCERS_ODE_ADVANCE_TO_TARGET_TIME_HPP_

#include "./impl/ode_advance_noop_observer.hpp"
#include "./impl/ode_advance_to_target_time.hpp"

namespace pressio{ namespace ode{

template<
  class StepperType,
  class StateType,
  class TimeType,
  class StepSizeSetterType,
  class ...Args
  >
void advance_to_target_time(StepperType & stepper,
			    StateType & ode_state,
			    const TimeType start_time,
			    const TimeType final_time,
			    StepSizeSetterType && dtManager,
			    Args&& ... args)
{

  // checks
  impl::static_assert_is_steppable_with(stepper, ode_state, start_time, std::forward<Args>(args)...);
  impl::static_assert_admissible_time_step_setter(dtManager, start_time);

  // run
  using observer_t = ::pressio::ode::impl::NoOpObserver<TimeType, StateType>;
  observer_t observer;
  impl::integrate_to_target_time_with_time_step_size_manager<false>
    (stepper, start_time, final_time, ode_state,
     observer, std::forward<StepSizeSetterType>(dtManager),
     std::forward<Args>(args)...);
}

template<
  class StepperType,
  class StateType,
  class TimeType,
  class StepSizeSetterType,
  class ObserverType,
  class ...Args
  >
::pressio::mpl::enable_if_t<
  ::pressio::ode::observer<ObserverType, TimeType, StateType>::value
  >
advance_to_target_time_and_observe(StepperType & stepper,
				   StateType & ode_state,
				   const TimeType	start_time,
				   const TimeType	final_time,
				   StepSizeSetterType	&& dtManager,
				   ObserverType	& observer,
				   Args&& ...args)
{

  // checks
  impl::static_assert_is_steppable_with(stepper, ode_state, start_time, std::forward<Args>(args)...);
  impl::static_assert_admissible_time_step_setter(dtManager, start_time);
  impl::static_assert_admissible_observer(observer, ode_state, start_time);

  // run
  impl::integrate_to_target_time_with_time_step_size_manager<false>
    (stepper, start_time, final_time, ode_state, observer,
     std::forward<StepSizeSetterType>(dtManager),
     std::forward<Args>(args)...);
}

template<
  class StepperType,
  class StateType,
  class TimeType,
  class StepSizeSetterType,
  class ObserverType,
  class ...Args
  >
void advance_to_target_time_with_time_step_recovery_and_observe(StepperType & stepper,
								StateType & ode_state,
								const TimeType start_time,
								const TimeType final_time,
								StepSizeSetterType && dtManager,
								ObserverType & observer,
								Args&& ... args)
{

  // checks
  impl::static_assert_is_steppable_with(stepper, ode_state, start_time, std::forward<Args>(args)...);
  impl::static_assert_admissible_time_step_setter(dtManager, start_time);
  impl::static_assert_admissible_observer(observer, ode_state, start_time);

  // run
  impl::integrate_to_target_time_with_time_step_size_manager<true>
    (stepper, start_time, final_time, ode_state,
     observer, std::forward<StepSizeSetterType>(dtManager),
     std::forward<Args>(args)...);
}

}}//end namespace pressio::ode
#endif  // ODE_ADVANCERS_ODE_ADVANCE_TO_TARGET_TIME_HPP_
