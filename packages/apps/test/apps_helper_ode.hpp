
#ifndef APPS_HELPER_ODE_HPP_
#define APPS_HELPER_ODE_HPP_

#include "ODE_ALL"
#include "SOLVERS_EXP"

namespace apps{
namespace test{

template <typename state_t,
	  typename residual_t,
	  typename scalar_t,
	  typename model_t,
	  typename time_t,
	  typename sizer_t,
	  typename aux_start_stepper_t,
	  typename residual_policy_t,
	  typename coll_t>
struct appsExpOdeHelper
{

  template<typename ... Args>
  static void explicitStandardEuler(state_t & y,
				    time_t ti,
				    time_t te,
				    time_t dt,
				    model_t & appObj,
				    coll_t & collectorObj,
				    Args && ... rest)
  {
    using stepper_t = ode::explicitEulerStepper<
      state_t, residual_t, scalar_t,
      model_t, scalar_t, sizer_t>;
    stepper_t stepperObj(appObj, std::forward<Args>(rest)... );

    ode::integrateNSteps(stepperObj, y, ti, dt, te/dt, collectorObj);    
  }//end 

  static void explicitStandardRK4(state_t & y,
				  time_t ti,
				  time_t te,
				  time_t dt,
				  model_t & appObj,
				  coll_t & collectorObj)
  {
    using stepper_t = ode::explicitRungeKutta4Stepper<
      state_t, residual_t, scalar_t,
      model_t, scalar_t, sizer_t>;
    stepper_t stepperObj(appObj);

    ode::integrateNSteps(stepperObj, y, ti, dt, te/dt, collectorObj);    
  }//end 

  
};//end class


//////////////////////////////////////////////
//////////////////////////////////////////////


template <typename state_t,
	  typename residual_t,
	  typename jacobian_t,
	  typename scalar_t,
	  typename model_t,
	  typename time_t,
	  typename sizer_t,
	  typename solver_policy_t,
	  typename aux_start_stepper_t,
	  typename residual_policy_t,
	  typename jacobian_policy_t,
	  typename coll_t>
struct appsImpOdeHelper
{
  // static void explicitStandardEuler(state_t & y,
  // 				    time_t ti,
  // 				    time_t te,
  // 				    time_t dt,
  // 				    model_t & appObj,
  // 				    coll_t & collectorObj)
  // {
  //   using stepper_t = ode::explicitEulerStepper<
  //     state_t, residual_t, scalar_t,
  //     model_t, scalar_t, sizer_t>;
  //   stepper_t stepperObj(appObj);
  //   ode::integrateNSteps(stepperObj, y, ti, dt, te/dt, collectorObj);    
  // }//end 
   
};//end class

  
}//end namespace test
}//end namespace apps
#endif 
