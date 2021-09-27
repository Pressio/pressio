/*
//@HEADER
// ************************************************************************
//
// ode_advance_n_steps.hpp
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

#ifndef ODE_ADVANCERS_IMPL_ODE_ADVANCE_N_STEPS_HPP_
#define ODE_ADVANCERS_IMPL_ODE_ADVANCE_N_STEPS_HPP_

#include "ode_advance_call_observer_dispatcher.hpp"
#include "ode_advance_printing_helpers.hpp"

namespace pressio{ namespace ode{ namespace impl{

template <
  class StepperType,
  class TimeType,
  class dt_setter,
  class StateType,
  class ObserverType,
  class StepCountType,
  class ... Args
  >
void advance_n_steps_with_dt_setter(StepperType & stepper,
				    const StepCountType & numSteps,
				    const TimeType & start_time,
				    StateType & odeStateInOut,
				    dt_setter && dtManager,
				    ObserverType & observer,
				    Args && ... args)
{
  using step_t = StepCountType;

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
  auto timer = Teuchos::TimeMonitor::getStackedTimer();
  timer->start("time loop");
#endif

  TimeType time = start_time;
  // pass initial condition to observer object
  call_observer(observer,
		 ::pressio::utils::Constants<step_t>::zero(),
		 time, odeStateInOut);

  TimeType dt = {};
  step_t step = 1;
  PRESSIOLOG_INFO("impl: advance_n_steps_with_dt_setter");
  for( ; step <= numSteps ; ++step)
    {
      // call the dt manager to set the dt to use
      dtManager(step, time, dt);
      print_step_time(step, time, dt);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer->start("time step");
#endif
      stepper(odeStateInOut, time, dt, step, std::forward<Args>(args)...);
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer->stop("time step");
#endif

      time += dt;
      call_observer(observer, step, time, odeStateInOut);
    }

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
  timer->stop("time loop");
#endif
}


template <
  class StepperType,
  class TimeType,
  class StateType,
  class ObserverType,
  class StepCountType,
  class ... Args
  >
void advance_n_steps_with_fixed_dt(StepperType & stepper,
				   const StepCountType & numSteps,
				   const TimeType & start_time,
				   const TimeType & step_size,
				   StateType & odeStateInOut,
				   ObserverType & observer,
				   Args && ... args)
{

  /* note that this function could be implemented
     in terms of the advance_with_dt_setter by simply
     passing a lambda that sets the dt to be same all the time.
     However, we don't do that because knowing we
     have a fixed dt, we can compute the time
     at each step by doing time = dt*step, whcih is
     better than just incrementing time for error accumulation.
  */

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
  auto timer = Teuchos::TimeMonitor::getStackedTimer();
  timer->start("time loop");
#endif

  using step_t = StepCountType;

  TimeType time = start_time;
  call_observer(observer,
		 ::pressio::utils::Constants<step_t>::zero(),
		 time, odeStateInOut);

  step_t step = 1;
  PRESSIOLOG_INFO("impl: advance_n_steps_with_fixed_dt");
  for( ; step <= numSteps ; ++step)
    {
      print_step_time(step, time, step_size);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer->start("time step");
#endif
      stepper(odeStateInOut, time, step_size, step, std::forward<Args>(args)...);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer->stop("time step");
#endif

      time = start_time + static_cast<TimeType>(step) * step_size;
      call_observer(observer, step, time, odeStateInOut);
    }
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
  timer->stop("time loop");
#endif
}


}}}//end namespace pressio::ode::impl
#endif  // ODE_ADVANCERS_IMPL_ODE_ADVANCE_N_STEPS_HPP_
