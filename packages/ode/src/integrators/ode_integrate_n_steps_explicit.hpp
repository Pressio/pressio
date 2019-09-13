/*
//@HEADER
// ************************************************************************
//
// ode_integrate_n_steps_explicit.hpp
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

#ifndef ODE_INTEGRATORS_INTEGRATE_N_STEPS_EXPLICIT_HPP_
#define ODE_INTEGRATORS_INTEGRATE_N_STEPS_EXPLICIT_HPP_

#include "ode_integrate_n_steps_impl.hpp"
#include "../meta/ode_is_legitimate_collector.hpp"

namespace pressio{ namespace ode{

template<
  typename stepper_type,
  typename state_type,
  typename time_type,
  typename integral_type,
  typename collector_type,
  typename std::enable_if<
    ode::meta::is_legitimate_collector<
      collector_type, integral_type, time_type, state_type
      >::value &&
    std::is_integral<integral_type>::value &&
    details::traits<stepper_type>::is_explicit
    >::type * = nullptr
  >
void integrateNSteps(stepper_type   & stepper,
		     state_type	    & yIn,
		     time_type	      start_time,
		     time_type	      dt,
		     integral_type    num_steps,
		     collector_type & collector)
{
  using empty_t = utils::impl::empty;
  using do_step_policy_t = impl::DoStepPolicy<empty_t, empty_t>;
  using advancer_t = impl::AdvancerPolicy<collector_type, do_step_policy_t>;
  advancer_t::execute(num_steps, start_time, dt, yIn, collector, stepper);
}

template<
  typename stepper_type,
  typename state_type,
  typename time_type,
  typename integral_type,
  typename std::enable_if<
    details::traits<stepper_type>::is_explicit
    >::type * = nullptr
  >
void integrateNSteps(stepper_type & stepper,
		     state_type	  & yIn,
		     time_type	    start_time,
		     time_type	    dt,
		     integral_type  num_steps){

  using empty_t = utils::impl::empty;
  using do_step_policy_t = impl::DoStepPolicy<empty_t, empty_t>;
  using advancer_t = impl::AdvancerPolicy<empty_t, do_step_policy_t>;
  advancer_t::execute(num_steps, start_time, dt, yIn, stepper);
}

}}//end namespace pressio::ode
#endif
