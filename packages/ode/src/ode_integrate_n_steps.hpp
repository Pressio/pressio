
#ifndef ODE_INTEGRATE_N_STEPS_HPP_
#define ODE_INTEGRATE_N_STEPS_HPP_

#include "ode_ConfigDefs.hpp"
#include "ode_forward_declarations.hpp"


template<typename stepper_type,
	 typename rhs_functor_type, 
	 typename state_type,
	 typename collector_functor_type
	 >
void integrate_n_steps<std::enable_if<!std::is_same<stepper_type,void>::value &&
				      !std::is_same<functor_type,void>::value 
				      >::type
		       >(stepper_type & stepper,
			 rhs_functor_type & rhs_functor,
			 state_type & start_state,
			 collector_functor_type & collector_fnctr,
			 time_type dt,
			 size_t num_of_steps)
{
  //typename time_type = typename core::defaultTypes::scalar_t;
  //typename odeint::unwrap_reference< Observer >::type &obs = observer;
  //typename odeint::unwrap_reference< Stepper >::type &st = stepper;

  using time_type = ode::details::time_type;
  time_type time = start_time;

  for( size_t step = 0; step < num_of_steps ; ++step )
  {
    collector_fnctr( start_state , time );
    stepper.do_step( rhs_functor, start_state, time, dt );
    time += + dt; //( step+1 ) * dt;
  }
  collector_fnctr( start_state , time );
}


#endif 
