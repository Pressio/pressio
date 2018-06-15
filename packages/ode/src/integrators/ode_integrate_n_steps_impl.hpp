
#ifndef ODE_INTEGRATE_N_STEPS_IMPL_HPP_
#define ODE_INTEGRATE_N_STEPS_IMPL_HPP_

#include "ode_ConfigDefs.hpp"
//#include "ode_forward_declarations.hpp"
#include "ode_void_collector.hpp"

namespace ode{
namespace impl{

template<typename stepper_type,
	 typename state_type,
	 typename time_type,
	 typename collector_type,
	 typename std::enable_if<
	   !std::is_void<stepper_type>::value &&
	   !std::is_same<collector_type,
			ode::voidCollector>::value
	   >::type * = nullptr
	 >
void integrateNStepsImpl(stepper_type & stepper,
			 state_type & stateIn,
			 time_type start_time,
			 time_type dt,
			 size_t num_of_steps,
			 collector_type & collector)
{
  time_type time = start_time;
  size_t step = 0;
  for( ; step < num_of_steps ; ++step)
  {
    // call collector/observer
    collector(step, time, stateIn);
    // do one step
    stepper.doStep(stateIn, time, dt);

    // advance time
    //mulitply (vs adding) benefits roundoff
    time = static_cast<double>(step) * dt;  
  }
  collector(step, time, stateIn);
}


//----------------------------------------------------

  
  
template<typename stepper_type,
	 typename state_type,
	 typename time_type,
	 typename std::enable_if<
	   !std::is_void<stepper_type>::value
	   >::type * = nullptr
	 >
void integrateNStepsImpl(stepper_type & stepper,
			 state_type & stateIn,
			 time_type start_time,
			 time_type dt,
			 size_t num_of_steps)
{
  time_type time = start_time;
  size_t step = 0;
  for( ; step < num_of_steps ; ++step)
  {
    // do one step
    stepper.doStep(stateIn, time, dt);
    // advance time
    //mulitply (vs adding) benefits roundoff
    time = static_cast<double>(step) * dt;  
  }
}
  

  
}//end namespace impl
}//end namespace ode

#endif 
