
#ifndef ODE_EXPLICIT_EULER_STEPPER_HPP_
#define ODE_EXPLICIT_EULER_STEPPER_HPP_

#include "./impl/ode_explicit_euler_stepper_impl.hpp"

namespace ode{

///////////////////////
// Standard policy 
/////////////////////// 
  
template<typename state_type,
	 typename residual_type,
	 typename scalar_type,
	 typename model_type,	
	 typename time_type,
	 typename sizer_type
	 >
class explicitEulerStepper<state_type,
			   residual_type,
			   scalar_type,
			   model_type,
			   time_type,
			   sizer_type,
			   void,
			   typename
			   std::enable_if<
			     !std::is_void<state_type>::value
			     >::type
			   >
  : public impl::explicitEulerStepperImpl<state_type,
					  residual_type,
					  scalar_type,
					  model_type,
					  time_type,
					  sizer_type,
					  ode::policy::residualStandardPolicy<
					    state_type, residual_type,
					    model_type, time_type, sizer_type>
					  >
{
public:
  using pol_t = ode::policy::residualStandardPolicy<
  state_type, residual_type, model_type, time_type, sizer_type>;

  using base_t = impl::explicitEulerStepperImpl<state_type,
						residual_type,
						scalar_type,
						model_type,
						time_type,
						sizer_type,
						pol_t>;
public:
  template < typename T = model_type,
	     typename... Args>
  explicitEulerStepper(T & model,
		       Args&&... rest)
    : base_t(model, policy_, std::forward<Args>(rest)...)
  {}
  explicitEulerStepper() = delete;
  ~explicitEulerStepper() = default;

private:
  pol_t policy_;

};//end class

  

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////



///////////////////////
// NON Standard policy 
/////////////////////// 
  
template<typename state_type,
	 typename residual_type,
	 typename scalar_type,
	 typename model_type,	
	 typename time_type,
	 typename sizer_type,
	 typename residual_policy_type
	 >
class explicitEulerStepper<state_type,
			   residual_type,
			   scalar_type,
			   model_type,
			   time_type,
			   sizer_type,
			   residual_policy_type,
			   typename
			   std::enable_if<
			     !std::is_void<residual_policy_type>::value
			     >::type
			   >
  : public impl::explicitEulerStepperImpl<state_type,
					  residual_type,
					  scalar_type,
					  model_type,
					  time_type,
					  sizer_type,
					  residual_policy_type>
{
public:
  using base_t = impl::explicitEulerStepperImpl<state_type,
						residual_type,
						scalar_type,
						model_type,
						time_type,
						sizer_type,
						residual_policy_type>;
public:
  template < typename T = model_type,
	     typename U = residual_policy_type,
	     typename... Args>
  explicitEulerStepper(T & model,
		       U & policy,
		       Args&&... rest)
    : base_t(model, policy, std::forward<Args>(rest)...)
  {}
  explicitEulerStepper() = delete;
  ~explicitEulerStepper() = default;

};//end class



}//end namespace
#endif 
