
#ifndef ODE_EXPLICIT_STEPPERS_EXPLICIT_EULER_STEPPER_USING_STD_POLICY_HPP_
#define ODE_EXPLICIT_STEPPERS_EXPLICIT_EULER_STEPPER_USING_STD_POLICY_HPP_

#include "./impl/ode_explicit_euler_stepper_impl.hpp"

namespace rompp{
namespace ode{

/////////////////////////////////////////
// State and resid are same type
/////////////////////////////////////////
 
// template<typename ode_state_type,
// 	 typename model_type>
// class ExplicitEulerStepper<ode_state_type, ode_state_type,
// 			   model_type, void,
// 			   typename
// 			   std::enable_if<
// 			     !std::is_void<ode_state_type>::value
// 			     >::type
// 			   >

struct empty{};


template<typename ode_state_type,
	 typename model_type>
class ExplicitEulerStepper<ExplicitSteppersEnum::Euler,
			   ode_state_type, 
			   model_type>
  : public std::conditional<ode::meta::isLegitimateExplicitStateType<ode_state_type>::value and
			 !ode::meta::is_legitimate_explicit_residual_policy<model_type>::value,
			    impl::ExplicitEulerStepperImpl<ode_state_type, ode_state_type,
			         typename core::details::traits<ode_state_type>::scalar_t,
				model_type, ode::policy::ExplicitResidualStandardPolicy<
					    ode_state_type, ode_state_type, model_type>>,
			    empty>::type
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


}//end namespace
}//end namespace rompp
#endif 
