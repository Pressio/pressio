
#ifndef ODE_INTEGRATORS_INTEGRATE_N_STEPS_EXPLICIT_HPP_
#define ODE_INTEGRATORS_INTEGRATE_N_STEPS_EXPLICIT_HPP_

#include "ode_integrate_n_steps_impl.hpp"
#include "../meta/ode_is_legitimate_collector.hpp"

namespace rompp{ namespace ode{

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
  using empty_t = core::impl::empty;
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

  using empty_t = core::impl::empty;
  using do_step_policy_t = impl::DoStepPolicy<empty_t, empty_t>;
  using advancer_t = impl::AdvancerPolicy<empty_t, do_step_policy_t>;
  advancer_t::execute(num_steps, start_time, dt, yIn, stepper);
}

}}//end namespace rompp::ode
#endif
