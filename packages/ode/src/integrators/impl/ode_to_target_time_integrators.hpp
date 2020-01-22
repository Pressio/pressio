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

#ifndef ODE_TO_TARGET_TIME_INTEGRATORS_HPP_
#define ODE_TO_TARGET_TIME_INTEGRATORS_HPP_

#include "../../ode_ConfigDefs.hpp"
#include "../../ode_fwd.hpp"

namespace pressio{ namespace ode{ namespace impl{

/*
 * time step size setter is passed by user
 * this is within the impl namespace, so should not be used outside
 */
template <typename DoStepPolicy_t>
struct IntegratorToTargetTimeWithTimeStepSizeSetter
{

  template <typename time_type, typename dt_setter, typename ... Args>
  static void execute(const time_type	& start_time,
		      const time_type	& final_time,
		      dt_setter		&& dtManager,
		      Args		&& ... args)
  {

    using step_t = ::pressio::ode::types::step_t;
    constexpr auto zero = ::pressio::utils::constants::zero<step_t>();

    if (final_time < start_time)
      throw std::runtime_error("You cannot call an integrator with a final time < start time.");

    if (final_time == start_time)
      return;

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("time loop");
#endif

    time_type time = start_time;
    time_type dt = {};

    step_t step	   = 1;
    ::pressio::utils::io::print_stdout("\nstarting time loop","\n");
    while (time < final_time)
    {
      // call the dt manager to set the dt to use for current step
      dtManager(step, time, dt);

#ifdef PRESSIO_ENABLE_DEBUG_PRINT
      auto fmt = utils::io::bg_grey() + utils::io::bold() + utils::io::red();
      auto reset = utils::io::reset();
      ::pressio::utils::io::print_stdout(fmt, "time step =", step, reset, "\n");
#endif

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer->start("time step");
#endif
      DoStepPolicy_t::execute(time, dt, step, std::forward<Args>(args)...);
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
      timer->stop("time step");
#endif

      time = start_time + static_cast<time_type>(step) * dt;
      step++;
    }

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("time loop");
#endif
  }//end ()
};


}}}//end namespace pressio::ode::impl
#endif
