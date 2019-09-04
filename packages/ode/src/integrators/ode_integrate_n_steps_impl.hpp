/*
//@HEADER
// ************************************************************************
//
// ode_integrate_n_steps_impl.hpp
//                     		      Pressio 
// Copyright 2019 National Technology & Engineering Solutions of Sandia,LLC 
//							      (NTESS)
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

#ifndef ODE_INTEGRATORS_INTEGRATE_N_STEPS_IMPL_HPP_
#define ODE_INTEGRATORS_INTEGRATE_N_STEPS_IMPL_HPP_

#include "../ode_ConfigDefs.hpp"
#include "../ode_fwd.hpp"

namespace pressio{ namespace ode{ namespace impl{

template< typename solver_type, typename guesser_cb_t>
struct DoStepPolicy{
  template <typename time_type,
	    typename integral_type,
	    typename state_type,
	    typename stepper_type>
  static void execute(time_type time,
		      time_type dt,
		      integral_type step,
		      state_type & yIn,
		      stepper_type & stepper,
		      solver_type & solver,
		      guesser_cb_t && guessCb)
  {
    stepper(yIn, time, dt, step, solver,
	    std::forward<guesser_cb_t>(guessCb));
  }
};

template<>
struct DoStepPolicy<utils::impl::empty, utils::impl::empty>{
  template <typename time_type,
	    typename integral_type,
	    typename state_type,
	    typename stepper_type>
  static void execute(time_type time,
		      time_type dt,
		      integral_type step,
		      state_type & yIn,
		      stepper_type & stepper)
  {
    stepper(yIn, time, dt, step);
  }
};


template<typename solver_type>
struct DoStepPolicy<solver_type, utils::impl::empty>{
  template <typename time_type,
	    typename integral_type,
	    typename state_type,
	    typename stepper_type>
  static void execute(time_type time,
		      time_type dt,
		      integral_type step,
		      state_type & yIn,
		      stepper_type & stepper,
		      solver_type & solver)
  {
    stepper(yIn, time, dt, step, solver);
  }
};
//-------------------------------------------------------

template <
  typename collector_type, typename int_type,
  typename time_type, typename state_type,
  typename enable = void
  >
struct CallCollectorDispatch;


template <
  typename collector_type, typename int_type,
  typename time_type, typename state_type
  >
struct CallCollectorDispatch<
  collector_type, int_type, time_type, state_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_wrapper<state_type>::value and
    ::pressio::ode::meta::collector_accepts_native_container
    <collector_type, int_type, time_type, state_type>::value
    >
  >
{
  static void execute(collector_type & collector,
		      const int_type & step,
		      const time_type & time,
		      const state_type & yIn){
    collector(step, time, *yIn.data());
  }
};


template <
  typename collector_type, typename int_type,
  typename time_type, typename state_type
  >
struct CallCollectorDispatch<
  collector_type, int_type, time_type, state_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_wrapper<state_type>::value and
    ::pressio::ode::meta::collector_accepts_pressio_container
    <collector_type, int_type, time_type, state_type>::value
    >
  >
{
  static void execute(collector_type & collector,
		      const int_type & step,
		      const time_type & time,
		      const state_type & yIn){
    collector(step, time, yIn);
  }
};



/*
 * A valid collector object is passed by user
 * to take snapshots
 */
template <typename collector_type, typename DoStepPolicy_t>
struct AdvancerPolicy{

  template <typename integral_type,
	    typename time_type,
	    typename state_type,
	    typename ... Args>
  static void execute(integral_type    num_steps,
		      time_type	   start_time,
		      time_type	   dt,
		      state_type &	   yIn,
		      collector_type & collector,
		      Args && ...	   args)
  {
#ifdef HAVE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("time loop");
#endif

    using collector_dispatch = CallCollectorDispatch<collector_type, integral_type,
						     time_type, state_type>;
    constexpr auto zero = ::pressio::utils::constants::zero<integral_type>();

    // time variable
    time_type time = start_time;
    // pass initial condition to collector object
    collector_dispatch::execute(collector, zero, time, yIn);
    //collector(0, time, yIn);

    integral_type step = 1;
    ::pressio::utils::io::print_stdout("\nstarting time loop","\n");
    for( ; step <= num_steps ; ++step)
    {
      #ifdef DEBUG_PRINT
      auto fmt = utils::io::bg_grey() + utils::io::bold() + utils::io::red();
      auto reset = utils::io::reset();
      ::pressio::utils::io::print_stdout(fmt, "time step =",
				      step, reset, "\n");
      #endif

#ifdef HAVE_TEUCHOS_TIMERS
      timer->start("time step");
#endif
      DoStepPolicy_t::execute(time, dt, step, yIn, std::forward<Args>(args)...);
#ifdef HAVE_TEUCHOS_TIMERS
      timer->stop("time step");
#endif

      time = start_time + static_cast<time_type>(step) * dt;
      collector_dispatch::execute(collector, step, time, yIn);
      //collector(step, time, yIn);
    }
#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("time loop");
#endif
  }//end ()
};



/*
 * No collector object is passed by user
 */
template <typename DoStepPolicy_t>
struct AdvancerPolicy<utils::impl::empty, DoStepPolicy_t>{

  template <typename integral_type,
	    typename time_type,
	    typename ... Args>
  static void execute(integral_type num_steps,
		      time_type	start_time,
		      time_type	dt,
		      Args && ... args)
  {
#ifdef HAVE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("time loop");
#endif

    time_type time = start_time;
    integral_type step = 1;

    ::pressio::utils::io::print_stdout("\nstarting time loop","\n");
    for( ; step <= num_steps ; ++step)
    {
      #ifdef DEBUG_PRINT
      auto fmt = utils::io::bg_grey() + utils::io::bold() + utils::io::red();
      auto reset = utils::io::reset();
      ::pressio::utils::io::print_stdout(fmt, "time step =",
				      step, reset, "\n");
      #endif

#ifdef HAVE_TEUCHOS_TIMERS
      timer->start("time step");
#endif
      DoStepPolicy_t::execute(time, dt, step, std::forward<Args>(args)...);
#ifdef HAVE_TEUCHOS_TIMERS
      timer->stop("time step");
#endif

      time = start_time + static_cast<time_type>(step) * dt;
    }

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("time loop");
#endif

  }//end ()
};


}}}//end namespace pressio::ode::impl
#endif
