
#ifndef ODE_INTEGRATE_N_STEPS_HPP_
#define ODE_INTEGRATE_N_STEPS_HPP_

#include "ode_ConfigDefs.hpp"
#include "ode_forward_declarations.hpp"


namespace ode{

  
template<typename stepper_type,
	 typename rhs_functor_type, 
	 typename state_type,
	 typename collector_functor_type = void
	 >

typename std::enable_if<!std::is_same<stepper_type,void>::value &&
			!std::is_same<rhs_functor_type,void>::value 
			>::type
integrate_n_steps(stepper_type & stepper,
		  rhs_functor_type & rhs_functor,
		  state_type & stateIn,
		  collector_functor_type & collector,
		  ode::details::time_type start_time,
		  ode::details::time_type dt,
		  size_t num_of_steps)
{
  //typename time_type = typename core::defaultTypes::scalar_t;
  //typename odeint::unwrap_reference< Observer >::type &obs = observer;
  //typename odeint::unwrap_reference< Stepper >::type &st = stepper;

  ode::details::time_type time = start_time;
  size_t step = 0;
  for( ; step < num_of_steps ; ++step)
  {
    // call collector/observer
    collector(step, time, stateIn);

    // do one step
    stepper.doStep( rhs_functor, stateIn, time, dt );

    // advance time
    time += + dt;
  }
  collector(step, time, stateIn);
}


}//end namespace

#endif 
