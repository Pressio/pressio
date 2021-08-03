/*
//@HEADER
// ************************************************************************
//
// ode_to_target_time_integrators.hpp
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

#ifndef ODE_INTEGRATORS_IMPL_ODE_TO_TARGET_TIME_INTEGRATORS_HPP_
#define ODE_INTEGRATORS_IMPL_ODE_TO_TARGET_TIME_INTEGRATORS_HPP_

#include "ode_advance_call_collector_dispatcher.hpp"
#include "ode_advance_printing_helpers.hpp"

namespace pressio{ namespace ode{ namespace impl{

template <bool b, typename TimeStepSizeManagerType, typename TimeType>
mpl::enable_if_t<b==true>
call_dt_manager(TimeStepSizeManagerType && dtManager,
	      const ::pressio::ode::step_count_type & step,
	      const TimeType & time,
	      TimeType & dt,
	      TimeType & minDt,
	      TimeType & dtRedFactor)
{
  dtManager(step, time, dt, minDt, dtRedFactor);
}

template <bool b, typename TimeStepSizeManagerType, typename TimeType>
mpl::enable_if_t<b==false>
call_dt_manager(TimeStepSizeManagerType && dtManager,
	      const ::pressio::ode::step_count_type & step,
	      const TimeType & time,
	      TimeType & dt,
	      TimeType & minDt,
	      TimeType & dtRedFactor)
{
  dtManager(step, time, dt);
}


template <
  bool enableTimeStepRecovery,
  typename StepperType,
  typename TimeType,
  typename StateType,
  typename CollectorType,
  typename TimeStepSizeManagerType,
  typename ... Args>
void
integrate_to_target_time_with_time_step_size_manager(StepperType & stepper,
					     const TimeType & start_time,
					     const TimeType & final_time,
					     StateType	& odeStateInOut,
					     CollectorType & collector,
					     TimeStepSizeManagerType	&& dtManager,
					     Args && ... args)
{

  using step_t = ::pressio::ode::step_count_type;
  constexpr auto zero = ::pressio::utils::Constants<step_t>::zero();

  if (final_time < start_time){
    throw std::runtime_error("You cannot call the advancer with final time < start time.");
  }

  if (final_time == start_time){
    return;
  }

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
  auto timer = Teuchos::TimeMonitor::getStackedTimer();
  timer->start("time loop");
#endif

  TimeType time  = start_time;
  TimeType dt    = pressio::utils::Constants<TimeType>::zero();
  TimeType minDt = pressio::utils::Constants<TimeType>::zero();
  TimeType dtReducFactor = ::pressio::utils::Constants<TimeType>::one();

  // pass initial condition to collector object
  call_collector(collector, zero, time, odeStateInOut);

  step_t step = 1;
  PRESSIOLOG_INFO("advance_to_target_timeWithDtCallback");
  constexpr auto eps = std::numeric_limits<TimeType>::epsilon();
  bool condition = true;
  while (condition)
    {
      PRESSIOLOG_DEBUG("callback dt manager");
      impl::call_dt_manager<enableTimeStepRecovery>(dtManager, step, time,
						  dt, minDt, dtReducFactor);


      if (dt <= static_cast<TimeType>(0)){
	throw std::runtime_error
	  ("The time step size cannot be <= 0.");
      }
      if (dt < minDt){
	throw std::runtime_error
	  ("The time step size cannot be smaller than the minimum value.");
      }

      if (enableTimeStepRecovery){
	if (minDt < static_cast<TimeType>(0)){
	  throw std::runtime_error
	    ("The minimum time step size cannot be smaller than zero.");
	}

	if (dtReducFactor <= ::pressio::utils::Constants<TimeType>::one()){
	  throw std::runtime_error
	    ("The time step size reduction factor must be > 1.");
	}
      }

      print_step_time(step, time, dt);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer->start("time step");
#endif

      if (enableTimeStepRecovery)
      {
	bool needStop = false;
	while(!needStop){
	  try
	  {
	    stepper.doStep(odeStateInOut, time, dt, step, std::forward<Args>(args)...);
	    needStop=true;
	  }
	  catch (::pressio::eh::TimeStepFailure const & e)
	  {
	    dt = dt/dtReducFactor;
	    if (dt<minDt){
	      throw std::runtime_error
		("Violation of minimum time step while trying to recover time step");
	    }

	    PRESSIOLOG_CRITICAL("time step={} failed, retrying with dt={}", step, dt);
	  }
	}
      }
      else
      {
	stepper.doStep(odeStateInOut, time, dt, step, std::forward<Args>(args)...);
      }

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer->stop("time step");
#endif

      time += dt;
      call_collector(collector, step, time, odeStateInOut);

      // use numeric limits to avoid tricky roundoff accumulation
      if ( std::abs(time - final_time) <= eps ) condition = false;

      // if we are over the final time, stop too
      if ( time > final_time ) condition = false;

      step++;
    }

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
  timer->stop("time loop");
#endif
}//end


}}}//end namespace pressio::ode::impl
#endif  // ODE_INTEGRATORS_IMPL_ODE_TO_TARGET_TIME_INTEGRATORS_HPP_
