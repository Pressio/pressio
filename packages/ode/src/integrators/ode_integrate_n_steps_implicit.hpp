/*
//@HEADER
// ************************************************************************
//
// ode_integrate_n_steps_implicit.hpp
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

#ifndef ODE_INTEGRATORS_INTEGRATE_N_STEPS_IMPLICIT_HPP_
#define ODE_INTEGRATORS_INTEGRATE_N_STEPS_IMPLICIT_HPP_

#include "../ode_ConfigDefs.hpp"
#include "./impl/ode_call_stepper_policy.hpp"
#include "./impl/ode_n_steps_integrators.hpp"
#include "../meta/ode_is_legitimate_collector.hpp"
#include "../meta/ode_is_legitimate_solver_for_implicit_stepper.hpp"

namespace pressio{ namespace ode{

// basic version
template<
  typename stepper_type,
  typename state_type,
  typename time_type,
  typename solver_type,
  typename std::enable_if<
    ::pressio::ode::details::traits<stepper_type>::is_implicit and
    !std::is_void<solver_type>::value and
    ::pressio::ode::meta::is_legitimate_solver_for_implicit_stepper<
      solver_type, stepper_type, state_type
      >::value
    >::type * = nullptr
  >
void integrateNSteps(stepper_type	 & stepper,
		     state_type		 & odeStateInOut,
		     const time_type	 start_time,
		     const time_type	 dt,
		     const types::step_t num_steps,
		     solver_type	 & solver)
{
  using empty_t = utils::impl::empty;
  using do_step_policy_t = impl::ImplicitDoStepBasic<solver_type>;
  using advancer_t = impl::IntegratorNStepsWithConstDt<do_step_policy_t>;
  advancer_t::execute(num_steps, start_time, dt, odeStateInOut, stepper, solver);
}


// with collector
template<
  typename stepper_type,
  typename state_type,
  typename time_type,
  typename collector_type,
  typename solver_type,
  typename std::enable_if<
    ::pressio::ode::details::traits<stepper_type>::is_implicit and
    ode::meta::is_legitimate_collector<
      collector_type, types::step_t, time_type, state_type
      >::value and
    ::pressio::ode::meta::is_legitimate_solver_for_implicit_stepper<
      solver_type, stepper_type, state_type
      >::value
    >::type * = nullptr
  >
void integrateNSteps(stepper_type	 & stepper,
		     state_type		 & odeStateInOut,
		     const time_type	 start_time,
		     const time_type	 dt,
		     const types::step_t num_steps,
		     collector_type	 & collector,
		     solver_type	 & solver)
{
  using empty_t = utils::impl::empty;
  using do_step_policy_t = impl::ImplicitDoStepBasic<solver_type>;
  using advancer_t = impl::IntegratorNStepsWithCollectorAndConstDt<collector_type, do_step_policy_t>;
  advancer_t::execute(num_steps, start_time, dt, odeStateInOut, collector, stepper, solver);
}


// with guesser
template<
  typename stepper_type,
  typename state_type,
  typename time_type,
  typename solver_type,
  typename guess_callback_t,
  typename std::enable_if<
    ::pressio::ode::details::traits<stepper_type>::is_implicit and
    ::pressio::ode::meta::is_legitimate_solver_for_implicit_stepper<
      solver_type, stepper_type, state_type
      >::value
    >::type * = nullptr
  >
void integrateNSteps(stepper_type	 & stepper,
		     state_type		 & odeStateInOut,
		     const time_type	 start_time,
		     const time_type	 dt,
		     const types::step_t num_steps,
		     solver_type	 & solver,
		     guess_callback_t    && guessCb)
{
  using empty_t = utils::impl::empty;
  using do_step_policy_t = impl::ImplicitDoStepWithGuesser<solver_type, guess_callback_t>;
  using advancer_t = impl::IntegratorNStepsWithConstDt<do_step_policy_t>;
  advancer_t::execute(num_steps, start_time, dt, odeStateInOut, stepper, solver,
		      std::forward<guess_callback_t>(guessCb));
}


// with guesser and collector
template<
  typename stepper_type,
  typename state_type,
  typename time_type,
  typename collector_type,
  typename solver_type,
  typename guess_callback_t,
  typename std::enable_if<
    ::pressio::ode::details::traits<stepper_type>::is_implicit and
    ::pressio::ode::meta::is_legitimate_collector<
      collector_type, types::step_t, time_type, state_type
      >::value and
    ::pressio::ode::meta::is_legitimate_solver_for_implicit_stepper<
      solver_type, stepper_type, state_type
      >::value
    >::type * = nullptr
  >
void integrateNSteps(stepper_type		& stepper,
		     state_type			& odeStateInOut,
		     const time_type	        start_time,
		     const time_type		dt,
		     const types::step_t	num_steps,
		     collector_type		& collector,
		     solver_type		& solver,
		     guess_callback_t		&& guessCb)
{
  using do_step_policy_t = impl::ImplicitDoStepWithGuesser<solver_type, guess_callback_t>;
  using advancer_t = impl::IntegratorNStepsWithCollectorAndConstDt<collector_type, do_step_policy_t>;
  advancer_t::execute(num_steps, start_time, dt, odeStateInOut, collector, stepper,
		      solver, std::forward<guess_callback_t>(guessCb));
}


// with guesser and collector passed in different order
template<
  typename stepper_type,
  typename state_type,
  typename time_type,
  typename collector_type,
  typename solver_type,
  typename guess_callback_t,
  typename std::enable_if<
    ::pressio::ode::meta::is_legitimate_collector<
      collector_type, types::step_t, time_type, state_type
      >::value and
    ::pressio::ode::meta::is_legitimate_solver_for_implicit_stepper<
      solver_type, stepper_type, state_type
      >::value
    >::type * = nullptr
  >
void integrateNSteps(stepper_type		& stepper,
		     state_type			& odeStateInOut,
		     const time_type	        start_time,
		     const time_type		dt,
		     const types::step_t	num_steps,
		     solver_type		& solver,
		     collector_type		& collector,
		     guess_callback_t		&& guessCb)
{
  integrateNSteps(stepper, odeStateInOut, start_time, dt, num_steps, collector, solver, guessCb);
}

}}//end namespace pressio::ode
#endif


// //---------------------------------------------------
// // arbitrary stepper with time step size scheduling
// //---------------------------------------------------
// template<
//   typename dt_setter_t,
//   typename stepper_type,
//   typename state_type,
//   typename time_type,
//   typename solver_type,
//   typename std::enable_if<
//     details::traits<stepper_type>::enum_id == ImplicitEnum::Arbitrary and
//     !std::is_void<dt_setter_t>::value
//     >::type * = nullptr
//   >
// void integrateNSteps(stepper_type   & stepper,
// 		     state_type	    & odeStateInOut,
// 		     const time_type	      start_time,
// 		     const ::pressio::ode::types::step_t    num_steps,
// 		     solver_type    & solver,
// 		     const dt_setter_t & dtManager)
// {
//   static constexpr bool with_dt = true;
//   using empty_t = utils::impl::empty;
//   using do_step_policy_t = impl::DoStepPolicy<solver_type, empty_t, with_dt>;
//   using advancer_t = impl::AdvancerNoCollectorForImplicitArbitraryStepper<do_step_policy_t>;
//   advancer_t::execute(num_steps, start_time, dtManager, odeStateInOut, stepper, solver);
// }
