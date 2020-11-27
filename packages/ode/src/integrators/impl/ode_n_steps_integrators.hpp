/*
//@HEADER
// ************************************************************************
//
// ode_n_steps_integrators.hpp
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

#ifndef ODE_INTEGRATORS_IMPL_ODE_N_STEPS_INTEGRATORS_HPP_
#define ODE_INTEGRATORS_IMPL_ODE_N_STEPS_INTEGRATORS_HPP_

#include "ode_call_collector_dispatcher.hpp"
#include "ode_integrators_printing_helpers.hpp"
#include "ode_dummy_collector.hpp"

namespace pressio{ namespace ode{ namespace impl{

struct IntegratorNStepsWithConstDt
{

  template <
    typename stepPolicy,
    typename time_type,
    typename state_type,
    typename collector_type,
    typename ... Args
    >
  static void execute(const ::pressio::ode::types::step_t & numSteps,
		      const time_type			  & start_time,
		      const time_type			  & dt,
		      state_type			  & odeStateInOut,
		      collector_type			  & collector,
		      Args				  && ... args)
  {
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("time loop");
#endif

    using step_t = ::pressio::ode::types::step_t;
    using collector_dispatch =
      CallCollectorDispatch<collector_type, time_type, state_type>;

    // time variable
    time_type time = start_time;
    // pass initial condition to collector object
    collector_dispatch::execute(collector,
				::pressio::utils::constants<step_t>::zero(),
				time, odeStateInOut);

    step_t step = 1;
    PRESSIOLOG_INFO("advanceNStepsWithConstDt");
    for( ; step <= numSteps ; ++step)
    {
      printStepTime(step, time, dt);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer->start("time step");
#endif
      stepPolicy::execute(time, dt, step, odeStateInOut,
			  std::forward<Args>(args)...);
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer->stop("time step");
#endif

      time = start_time + static_cast<time_type>(step) * dt;
      collector_dispatch::execute(collector, step, time, odeStateInOut);
    }
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("time loop");
#endif
  }//end ()
};


/*
 * time step size setter is passed
 */
struct IntegratorNStepsWithTimeStepSizeSetter
{

  template <
    typename stepPolicy,
    typename time_type,
    typename dt_setter,
    typename state_type,
    typename collector_type,
    typename ... Args
    >
  static void execute(const ::pressio::ode::types::step_t & numSteps,
		      const time_type			  & start_time,
		      state_type			  & odeStateInOut,
		      dt_setter				  && dtManager,
		      collector_type			  & collector,
		      Args				  && ... args)
  {
    using step_t = ::pressio::ode::types::step_t;
    using collector_dispatch =
      CallCollectorDispatch<collector_type, time_type, state_type>;

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("time loop");
#endif

    time_type time = start_time;
    // pass initial condition to collector object
    collector_dispatch::execute(collector,
				::pressio::utils::constants<step_t>::zero(),
				time, odeStateInOut);

    time_type dt = {};
    step_t step	   = 1;
    PRESSIOLOG_INFO("advanceNStepsWithDtCallBack");
    for( ; step <= numSteps ; ++step)
    {
      // call the dt manager to set the dt to use at the beginning
      dtManager(step, time, dt);
      printStepTime(step, time, dt);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer->start("time step");
#endif
      stepPolicy::execute(time, dt, step,
			  odeStateInOut, std::forward<Args>(args)...);
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer->stop("time step");
#endif

      time += dt;
      collector_dispatch::execute(collector, step, time, odeStateInOut);
    }

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("time loop");
#endif
  }
};

}}}//end namespace pressio::ode::impl
#endif  // ODE_INTEGRATORS_IMPL_ODE_N_STEPS_INTEGRATORS_HPP_
