
#ifndef ODE_STEPPERS_EXPLICIT_STEPPERS_EXPLICIT_EULER_STEPPER_HPP_
#define ODE_STEPPERS_EXPLICIT_STEPPERS_EXPLICIT_EULER_STEPPER_HPP_

#include "./impl/ode_explicit_euler_stepper_impl.hpp"

namespace rompp{
namespace ode{

/////////////////////////////////////////
//  (a) State and resid are same type
//  (b) Standard policy 
/////////////////////////////////////////
 
template<typename ode_state_type,
	 typename model_type>
class ExplicitEulerStepper<ode_state_type, ode_state_type,
			   model_type, void,
			   typename
			   std::enable_if<
			     !std::is_void<ode_state_type>::value
			     >::type
			   >
  : public impl::ExplicitEulerStepperImpl<
	     ode_state_type,
             ode_state_type,
	     typename core::details::traits<ode_state_type>::scalar_t,
	     model_type,
	     ode::policy::ExplicitResidualStandardPolicy<
	       ode_state_type, ode_state_type, model_type>
	     >
{

public:
  using pol_t = ode::policy::ExplicitResidualStandardPolicy<
     ode_state_type, ode_state_type, model_type>;

  using base_t = impl::ExplicitEulerStepperImpl<
     ode_state_type, ode_state_type,
     typename core::details::traits<ode_state_type>::scalar_t,
     model_type, pol_t>;

public:
  
  template < typename T1 = model_type,
	     typename T2 = ode_state_type,
	     typename... Args>
  ExplicitEulerStepper(T1 & model,
		       T2 const & y0,
		       T2 const & r0,
		       Args&&... rest)
    : base_t(model, policy_, y0, r0, std::forward<Args>(rest)...)
  {}

  ExplicitEulerStepper() = delete;
  ~ExplicitEulerStepper() = default;

private:
  pol_t policy_;

};//end class



/////////////////////////////////////////
//  (a) State and resid are same type
//  (b) NON Standard policy 
/////////////////////////////////////////
  
template<typename ode_state_type,
	 typename model_type,
	 typename residual_policy_type>
class ExplicitEulerStepper<ode_state_type,
			   ode_state_type,
			   model_type,
			   residual_policy_type,
			   typename
			   std::enable_if<
			     !std::is_void<residual_policy_type>::value 
			     >::type
			   >
  : public impl::ExplicitEulerStepperImpl<
	     ode_state_type,
	     ode_state_type,
	     typename core::details::traits<ode_state_type>::scalar_t,
	     model_type,
	     residual_policy_type>
{

public:
  using base_t = impl::ExplicitEulerStepperImpl<
     ode_state_type, ode_state_type,
     typename core::details::traits<ode_state_type>::scalar_t,
     model_type, residual_policy_type>;

public:
  template < typename T1 = model_type,
	     typename T2 = residual_policy_type,
	     typename T3 = ode_state_type,
	     typename... Args>
  ExplicitEulerStepper(T1 & model,
		       T2 & policy,
		       T3 const & y0,
		       T3 const & r0,
		       Args&&... rest)
    : base_t(model, policy, y0, r0,
	     std::forward<Args>(rest)...){}

  ExplicitEulerStepper() = delete;
  ~ExplicitEulerStepper() = default;

};//end class


}//end namespace
}//end namespace rompp
#endif 
