
#ifndef ODE_STEPPERS_EXPLICIT_STEPPERS_EXPLICIT_EULER_STEPPER_HPP_
#define ODE_STEPPERS_EXPLICIT_STEPPERS_EXPLICIT_EULER_STEPPER_HPP_

#include "./impl/ode_explicit_euler_stepper_impl.hpp"

namespace ode{

///////////////////////
// Standard policy 
/////////////////////// 
  
template<typename state_type,
	 typename space_residual_type,
	 typename scalar_type,
	 typename model_type,	
	 typename sizer_type
	 >
class ExplicitEulerStepper<state_type,
			   space_residual_type,
			   scalar_type,
			   model_type,
			   sizer_type,
			   void,
			   typename
			   std::enable_if<
			     !std::is_void<state_type>::value
			     >::type
			   >
  : public impl::ExplicitEulerStepperImpl<
		   state_type,
		   space_residual_type,
		   scalar_type,
		   model_type,
		   sizer_type,
		   ode::policy::explicit_residual_standard_policy<
		     state_type, space_residual_type,
		     model_type, scalar_type, sizer_type>
		   >
{

public:
  using pol_t = ode::policy::explicit_residual_standard_policy<
  state_type, space_residual_type, model_type, scalar_type, sizer_type>;

  using base_t = impl::ExplicitEulerStepperImpl<state_type,
						space_residual_type,
						scalar_type,
						model_type,
						sizer_type,
						pol_t>;
public:
  template < typename T1 = model_type,
	     typename T2 = state_type,
	     typename T3 = space_residual_type,
	     typename... Args>
  ExplicitEulerStepper(T1 & model, T2 const & y0,
		       T3 const & r0, Args&&... rest)
    : base_t(model, policy_, y0, r0, std::forward<Args>(rest)...)
  {}

  ExplicitEulerStepper() = delete;
  ~ExplicitEulerStepper() = default;

private:
  pol_t policy_;

};//end class

  

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////


///////////////////////
// NON Standard policy 
/////////////////////// 
  
template<typename state_type,
	 typename space_residual_type,
	 typename scalar_type,
	 typename model_type,	
	 typename sizer_type,
	 typename residual_policy_type
	 >
class ExplicitEulerStepper<state_type,
			   space_residual_type,
			   scalar_type,
			   model_type,
			   sizer_type,
			   residual_policy_type,
			   typename
			   std::enable_if<
			     !std::is_void<residual_policy_type>::value
			     >::type
			   >
  : public impl::ExplicitEulerStepperImpl<state_type,
					  space_residual_type,
					  scalar_type,
					  model_type,
					  sizer_type,
					  residual_policy_type>
{
public:
  using base_t = impl::ExplicitEulerStepperImpl<state_type,
						space_residual_type,
						scalar_type,
						model_type,
						sizer_type,
						residual_policy_type>;
public:
  template < typename T1 = model_type,
	     typename T2 = residual_policy_type,
	     typename T3 = state_type,
	     typename T4 = space_residual_type,
	     typename... Args>
  ExplicitEulerStepper(T1 & model,
		       T2 & policy,
		       T3 const & y0,
		       T4 const & r0,
		       Args&&... rest)
    : base_t(model, policy, y0, r0, std::forward<Args>(rest)...)
  {}

  ExplicitEulerStepper() = delete;
  ~ExplicitEulerStepper() = default;

};//end class


}//end namespace
#endif 
