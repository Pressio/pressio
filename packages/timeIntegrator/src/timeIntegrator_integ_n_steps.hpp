
#ifndef TIMEINTEGRATOR_INTEG_N_STEPS_HPP
#define TIMEINTEGRATOR_INTEG_N_STEPS_HPP

#include "timeIntegrator_ConfigDefs.hpp"

template<
  class stepper_type,
  class system_type, 
  class state_type,
  class time_type = typename timeIntegrator::details::defaultTypes::scalar_t,
  class observer_type = void
  >
time_type integrate_n_steps(stepper_type & stepper,
			    system_type & system,
			    state_type & start_state,
			    time_type start_time,
			    time_type dt,
			    size_t num_of_steps)
{
  //typename odeint::unwrap_reference< Observer >::type &obs = observer;
  //typename odeint::unwrap_reference< Stepper >::type &st = stepper;

  time_type time = start_time;

  for( size_t step = 0; step < num_of_steps ; ++step ){
    //obs( start_state , time );
    stepper.do_step( system, start_state, time, dt );
    time = start_time + static_cast<time_type>( step+1 ) * dt;
  }
  //  obs( start_state , time );
  return time;
}


#endif 
