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

#include "ode_call_collector_dispatcher.hpp"
#include "ode_integrators_printing_helpers.hpp"

namespace pressio{ namespace ode{ namespace impl{

template <bool b, typename dt_manager, typename time_type>
mpl::enable_if_t<b==true>
callDtManager(dt_manager && dtManager,
	      const ::pressio::ode::types::step_t & step,
	      const time_type & time,
	      time_type & dt,
	      time_type & minDt,
	      time_type & dtRedFactor)
{
  dtManager(step, time, dt, minDt, dtRedFactor);
}

template <bool b, typename dt_manager, typename time_type>
mpl::enable_if_t<b==false>
callDtManager(dt_manager && dtManager,
	      const ::pressio::ode::types::step_t & step,
	      const time_type & time,
	      time_type & dt,
	      time_type & minDt,
	      time_type & dtRedFactor)
{
  dtManager(step, time, dt);
}


template <
  bool enableTimeStepRecovery,
  typename stepPolicy,
  typename time_type,
  typename collector_t,
  typename dt_manager,
  typename state_type,
  typename ... Args>
void
integrateToTargetTimeWithTimeStepSizeManager(const time_type	& start_time,
					     const time_type	& final_time,
					     collector_t	& collector,
					     dt_manager	&& dtManager,
					     state_type	& odeStateInOut,
					     Args		&& ... args)
{

  using step_t = ::pressio::ode::types::step_t;
  using collector_dispatch = CallCollectorDispatch<collector_t, time_type, state_type>;
  constexpr auto zero = ::pressio::utils::constants<step_t>::zero();

  if (final_time < start_time)
    throw std::runtime_error("You cannot call the advancer with final time < start time.");

  if (final_time == start_time)
    return;

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
  auto timer = Teuchos::TimeMonitor::getStackedTimer();
  timer->start("time loop");
#endif

  time_type time  = start_time;
  time_type dt    = pressio::utils::constants<time_type>::zero();
  time_type minDt = pressio::utils::constants<time_type>::zero();
  time_type dtReducFactor = ::pressio::utils::constants<time_type>::one();

  // pass initial condition to collector object
  collector_dispatch::execute(collector, zero, time, odeStateInOut);

  step_t step = 1;
  PRESSIOLOG_INFO("advanceToTargetTimeWithDtCallback");
  constexpr auto eps = std::numeric_limits<time_type>::epsilon();
  bool condition = true;
  while (condition)
    {
      PRESSIOLOG_DEBUG("callback dt manager");
      impl::callDtManager<enableTimeStepRecovery>(dtManager, step, time,
						  dt, minDt, dtReducFactor);


      if (dt <= static_cast<time_type>(0)){
	throw std::runtime_error
	  ("The time step size cannot be <= 0.");
      }
      if (dt < minDt){
	throw std::runtime_error
	  ("The time step size cannot be smaller than the minimum value.");
      }

      if (enableTimeStepRecovery){
	if (minDt < static_cast<time_type>(0)){
	  throw std::runtime_error
	    ("The minimum time step size cannot be smaller than zero.");
	}

	if (dtReducFactor <= ::pressio::utils::constants<time_type>::one()){
	  throw std::runtime_error
	    ("The time step size reduction factor must be > 1.");
	}
      }

      printStepTime(step, time, dt);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer->start("time step");
#endif

      if (enableTimeStepRecovery)
      {
	bool needStop = false;
	while(!needStop){
	  try
	  {
	    stepPolicy::execute(time, dt, step, odeStateInOut,
				std::forward<Args>(args)...);
	    needStop=true;
	  }
	  catch (::pressio::eh::time_step_failure const & e)
	  {
	    dt = dt/dtReducFactor;
	    if (dt<minDt){
	      throw std::runtime_error
		("Violation of minimum time step while trying to recover time step");
	    }

	    PRESSIOLOG_CRITICAL("time step={} failed, retrying with dt={}", step, dt);
	    // auto fmt = ::pressio::utils::io::red();
	    // auto rst = ::pressio::utils::io::reset();
	    // ::pressio::utils::io::print_stdout
	    // 	(std::setw(1), fmt,
	    // 	 "time step=", step, "failed, retrying with dt=", dt,
	    // 	 rst, "\n");
	  }
	}
      }
      else
      {
	stepPolicy::execute(time, dt, step, odeStateInOut,
			    std::forward<Args>(args)...);
      }

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer->stop("time step");
#endif

      time += dt;
      collector_dispatch::execute(collector, step, time, odeStateInOut);

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
