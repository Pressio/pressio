
#ifndef ODE_INTEGRATORS_INTEGRATE_N_STEPS_IMPLICIT_HPP_
#define ODE_INTEGRATORS_INTEGRATE_N_STEPS_IMPLICIT_HPP_

#include "ode_integrate_n_steps_impl.hpp"
#include "../meta/ode_basic_meta.hpp"
#include "../meta/ode_is_legitimate_collector.hpp"

namespace rompp{ namespace ode{

template<typename stepper_type,   typename state_type,
	 typename time_type,	  typename integral_type,
	 typename collector_type, typename solver_type,
	 typename std::enable_if<
	   ode::meta::is_legitimate_collector<
	     collector_type, integral_type,
	     time_type, state_type>::value &&
	   std::is_integral<integral_type>::value &&
	   details::traits<stepper_type>::is_implicit
	   >::type * = nullptr>
void integrateNSteps(stepper_type   & stepper,
		     state_type	    & yIn,
		     time_type	      start_time,
		     time_type	      dt,
		     integral_type    num_steps,
		     collector_type & collector,
		     solver_type    & solver){
  using do_step_policy_t = impl::DoStepPolicy<solver_type,
					      core::impl::empty>;
  using advancer_t = impl::AdvancerPolicy<collector_type,
					  do_step_policy_t>;
  advancer_t advancer;
  advancer(num_steps, start_time, dt, yIn, collector, stepper, solver);
}


template<typename stepper_type,	 typename state_type,
	 typename time_type,	 typename integral_type,
	 typename solver_type,
	 typename std::enable_if<
	    details::traits<stepper_type>::is_implicit
	   >::type * = nullptr>
void integrateNSteps(stepper_type   & stepper,
		     state_type	    & yIn,
		     time_type	      start_time,
		     time_type	      dt,
		     integral_type    num_steps,
		     solver_type    & solver){

  using do_step_policy_t = impl::DoStepPolicy<solver_type,
					      core::impl::empty>;
  using advancer_t = impl::AdvancerPolicy<core::impl::empty,
					  do_step_policy_t>;
  advancer_t advancer;
  advancer(num_steps, start_time, dt, yIn, stepper, solver);
}


template<typename stepper_type, typename state_type,
	 typename time_type,	typename integral_type,
	 typename solver_type,  typename guess_callback_t,
	 typename std::enable_if<
	    details::traits<stepper_type>::is_implicit
	   >::type * = nullptr>
void integrateNSteps(stepper_type   & stepper,
		     state_type	    & yIn,
		     time_type	      start_time,
		     time_type	      dt,
		     integral_type    num_steps,
		     solver_type    & solver,
		     guess_callback_t   && guessCb){

  using do_step_policy_t = impl::DoStepPolicy<solver_type, guess_callback_t>;
  using advancer_t = impl::AdvancerPolicy<core::impl::empty,
					  do_step_policy_t>;
  advancer_t advancer;
  advancer(num_steps, start_time, dt, yIn, stepper, solver,
	   std::forward<guess_callback_t>(guessCb));
}


template<typename stepper_type,   typename state_type,
	 typename time_type,	  typename integral_type,
	 typename collector_type, typename solver_type,
	 typename guess_callback_t,
	 typename std::enable_if<
	   ode::meta::is_legitimate_collector<
	     collector_type, integral_type,
	     time_type, state_type>::value &&
	   std::is_integral<integral_type>::value &&
	   details::traits<stepper_type>::is_implicit
	   >::type * = nullptr>
void integrateNSteps(stepper_type   & stepper,
		     state_type	    & yIn,
		     time_type	      start_time,
		     time_type	      dt,
		     integral_type    num_steps,
		     collector_type & collector,
		     solver_type    & solver,
		     guess_callback_t && guessCb){
  using do_step_policy_t = impl::DoStepPolicy<solver_type, guess_callback_t>;
  using advancer_t = impl::AdvancerPolicy<collector_type,
					  do_step_policy_t>;
  advancer_t advancer;
  advancer(num_steps, start_time, dt, yIn, collector, stepper,
	   solver, std::forward<guess_callback_t>(guessCb));
}

}}//end namespace rompp::ode
#endif





// namespace meta{

// template<typename stepper_type,  typename state_type,
// 	 typename time_type,     typename integral_type,
// 	 typename collector_type = void,
// 	 typename enable = void>
// struct are_legitimate_types_for_nsteps_integration : std::false_type{};

// template<typename stepper_type,  typename state_type,
// 	 typename time_type,	 typename integral_type,
// 	 typename collector_type,
// 	 typename enable = void>
// struct are_legitimate_types_for_nsteps_integration<
//   stepper_type, state_type, time_type, integral_type,
//   core::meta::enable_if_t<
//     ode::meta::is_legitimate_collector<collector_type, integral_type,
// 				       time_type, state_type>::value &&
//     std::is_integral<integral_type>::value &&
//     >
//   > : std::true_type{};

// }// end namespace meta
