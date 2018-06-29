
#ifndef ODE_INTEGRATE_N_STEPS_HPP_
#define ODE_INTEGRATE_N_STEPS_HPP_

#include "ode_ConfigDefs.hpp"
#include "../meta/ode_meta.hpp"

namespace ode{
  
template<typename stepper_type,
	 typename state_type,
	 typename collector_type,
	 typename time_type,
	 typename integral_type,
	 typename std::enable_if<
	   !std::is_void<stepper_type>::value &&
	   ode::meta::isLegitimateCollector<collector_type,
					    integral_type,
					    time_type,
					    state_type>::value
	   >::type * = nullptr
	 >
void integrateNSteps(stepper_type & stepper,
		     state_type & stateIn,
		     collector_type & collector, 
		     time_type start_time,
		     time_type dt,
		     integral_type num_steps)
{
  time_type time = start_time;
  integral_type step = 0;

  // time loop
  for( ; step < num_steps ; ++step)
  {
    // call collector/observer 
    collector(step, time, stateIn);
    // do one step
    stepper.doStep(stateIn, time, dt);
    // advance time: mulitply (vs adding) benefits roundoff
    time = start_time + static_cast<time_type>(step) * dt;  
  }
  collector(step, time, stateIn);
}
//----------------------------------------------------------------


template<typename stepper_type,
	 typename state_type,
	 typename time_type,
	 typename integral_type,
	 typename std::enable_if<
	   !std::is_void<stepper_type>::value
	   >::type * = nullptr
	 >
void integrateNSteps(stepper_type & stepper,
		     state_type & stateIn,
		     time_type start_time,
		     time_type dt,
		     integral_type num_steps)
{
  time_type time = start_time;
  integral_type step = 0;
  for( ; step < num_steps ; ++step)
  {
    // do one step
    stepper.doStep(stateIn, time, dt);
    // advance time: mulitply (vs adding) benefits roundoff
    time = start_time + static_cast<time_type>(step) * dt;  
  }
}

}//end namespace
#endif 
