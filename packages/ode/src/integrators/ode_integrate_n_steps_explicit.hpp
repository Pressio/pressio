
#ifndef ODE_INTEGRATORS_INTEGRATE_N_STEPS_EXPLICIT_HPP_
#define ODE_INTEGRATORS_INTEGRATE_N_STEPS_EXPLICIT_HPP_

#include "ode_integrate_n_steps_impl.hpp"
#include "../meta/ode_basic_meta.hpp"
#include "../meta/ode_is_legitimate_collector.hpp"

namespace rompp{ namespace ode{

// template<typename stepper_type, typename state_type,
// 	 typename time_type, typename integral_type,
// 	 typename collector_type,
// 	 typename std::enable_if<
// 	   ode::meta::is_legitimate_collector<
// 	     collector_type, integral_type, time_type, state_type
// 	     >::value &&
// 	   std::is_integral<integral_type>::value &&
// 	   details::traits<stepper_type>::is_explicit
// 	   >::type * = nullptr>
// void integrateNSteps(stepper_type   & stepper,
// 		     state_type	    & yIn,
// 		     time_type	      start_time,
// 		     time_type	      dt,
// 		     integral_type    num_steps,
// 		     collector_type & collector)
// {
//   impl::integrateNSteps<stepper_type, state_type, time_type,
// 			integral_type, false,
// 			collector_type>(stepper, yIn, start_time,
// 					dt, num_steps, &collector);
// }

template<typename stepper_type, typename state_type,
	 typename time_type, typename integral_type,
	 typename std::enable_if<
	   details::traits<stepper_type>::is_explicit
	   >::type * = nullptr>
void integrateNSteps(stepper_type & stepper,
		     state_type	  & yIn,
		     time_type	    start_time,
		     time_type	    dt,
		     integral_type  num_steps){

  using do_step_policy_t = impl::DoStepPolicy<core::impl::empty,
					      core::impl::empty>;
  using advancer_t = impl::AdvancerPolicy<core::impl::empty,
					       do_step_policy_t>;
  advancer_t advancer;
  advancer(num_steps, start_time, dt, yIn, stepper);
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
//   ::rompp::mpl::enable_if_t<
//     ode::meta::is_legitimate_collector<collector_type, integral_type,
// 				       time_type, state_type>::value &&
//     std::is_integral<integral_type>::value &&
//     >
//   > : std::true_type{};

// }// end namespace meta
