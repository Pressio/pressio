
#ifndef ODE_INTEGRATE_N_STEPS_IMPL_HPP_
#define ODE_INTEGRATE_N_STEPS_IMPL_HPP_

#include "ode_ConfigDefs.hpp"

namespace ode{
namespace impl{
  
// template<typename stepper_type,
// 	 typename state_type,
// 	 typename time_type,
// 	 typename integral_type,
// 	 typename collector_type>
// void integrateNStepsImpl(stepper_type & stepper,
// 			 state_type & y,
// 			 time_type & start_time,
// 			 time_type dt,
// 			 integral_type num_steps,
// 			 collector_type & collector)
// {
//   time_type time = start_time;
//   integral_type step = 1;

//   // time loop
//   for( ; step <= num_steps ; ++step)
//   {
//     // call collector/observer 
//     collector(step, time, y);
//     // do one step
//     stepper.doStep(y, time, dt, step);
//     // advance time: mulitply (vs adding) benefits roundoff
//     time = start_time + static_cast<time_type>(step) * dt;  
//   }
//   collector(step, time, y);
// }


}//end namespace impl
}//end namespace ode
#endif 
