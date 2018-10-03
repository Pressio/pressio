
#ifndef ODE_EXPLICIT_STEPPERS_EXPLICIT_EULER_STEPPER_USING_ARBIT_POLICY_HPP_
#define ODE_EXPLICIT_STEPPERS_EXPLICIT_EULER_STEPPER_USING_ARBIT_POLICY_HPP_

#include "./impl/ode_explicit_euler_stepper_impl.hpp"

namespace rompp{
namespace ode{

/////////////////////////////////////////
// State and resid are same type
/////////////////////////////////////////
  
template<typename ode_state_type,
	 typename model_type,
	 typename residual_policy_type>
class ExplicitEulerStepper<ExplicitSteppersEnum::Euler,
			   ode_state_type,
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
