
#ifndef ODE_STEPPERS_IMPLICIT_STEPPERS_IMPLICIT_EULER_STEPPER_HPP_
#define ODE_STEPPERS_IMPLICIT_STEPPERS_IMPLICIT_EULER_STEPPER_HPP_

#include "./impl/ode_implicit_euler_stepper_impl.hpp"
#include "../../policies/standard/ode_implicit_euler_residual_standard_policy.hpp"
#include "../../policies/standard/ode_implicit_euler_jacobian_standard_policy.hpp"

namespace ode{

//*********************************************************
// residual policy = standard 
// jacobian policy = standard 
//*********************************************************
template<typename state_type,
	 typename residual_type,
	 typename jacobian_type, 
	 typename model_type>
class ImplicitEulerStepper<state_type, residual_type,
			   jacobian_type, 
			   model_type, 
			   void, void,
			   typename std::enable_if<
			     !std::is_void<state_type>::value
			     >::type>
  : public impl::ImplicitEulerStepperImpl<
       state_type, residual_type, jacobian_type,
       typename core::details::traits<state_type>::scalar_t,
       model_type, 
       ode::policy::ImplicitEulerResidualStandardPolicy<
	 state_type, residual_type, model_type>,
       ode::policy::ImplicitEulerJacobianStandardPolicy<
	 state_type, jacobian_type, model_type>
  >
{
private:
  using scalar_type = typename core::details::traits<state_type>::scalar_t;

public:
  using res_pol_t = ode::policy::ImplicitEulerResidualStandardPolicy<
  state_type, residual_type, model_type>;
  
  using jac_pol_t = ode::policy::ImplicitEulerJacobianStandardPolicy<
    state_type, jacobian_type, model_type>;

  using base_t = impl::ImplicitEulerStepperImpl<state_type,
						residual_type,
  						jacobian_type,
						scalar_type,
  						model_type,
  						res_pol_t,
						jac_pol_t>;
public:
  template < typename T1 = model_type,
	     typename T4 = state_type,
	     typename... Args>
  ImplicitEulerStepper(T1 & model,
		       T4 const & y0,
		       Args&&... rest)
    : base_t(model, res_policy_obj_, jac_policy_obj_,
	     y0, std::forward<Args>(rest)...)
  {}

  
  ImplicitEulerStepper() = delete;
  ~ImplicitEulerStepper() = default;

private:
  res_pol_t res_policy_obj_;
  jac_pol_t jac_policy_obj_;
  
}; //end class
  

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
  

  
//*********************************************************
// residual policy = arbitrary, non void
// jacobian policy = arbitrary, non void
//*********************************************************
template<typename state_type,
	 typename residual_type,
	 typename jacobian_type,
	 typename model_type,
	 typename residual_policy_type,
	 typename jacobian_policy_type
	 >
class ImplicitEulerStepper<state_type,
			   residual_type,
			   jacobian_type,
			   model_type,
			   residual_policy_type,
			   jacobian_policy_type,
			   typename
			   std::enable_if<
			     !std::is_void<residual_policy_type>::value &&
			     !std::is_void<jacobian_policy_type>::value
			     >::type
			   >
  : public impl::ImplicitEulerStepperImpl<state_type,
					  residual_type,
					  jacobian_type,
    typename core::details::traits<state_type>::scalar_t,
					  model_type,
					  residual_policy_type,
					  jacobian_policy_type>
{

private:
  using scalar_type = typename core::details::traits<state_type>::scalar_t;

public:
  using base_t = impl::ImplicitEulerStepperImpl<state_type,
						residual_type,
						jacobian_type,
						scalar_type,
						model_type,
						residual_policy_type,
						jacobian_policy_type>;
public:
  template < typename T1 = model_type,
	     typename T2 = residual_policy_type,
	     typename T3 = jacobian_policy_type,
	     typename T4 = state_type,
	     typename... Args>
  ImplicitEulerStepper(T1 & model,
		       T2 & res_policy_obj,
		       T3 & jac_policy_obj,
		       T4 const & y0,
		       Args&&... rest)
    : base_t(model, res_policy_obj, jac_policy_obj,
	     y0, std::forward<Args>(rest)...)
  {}
  
  ImplicitEulerStepper() = delete;
  ~ImplicitEulerStepper() = default;

}; //end class


}//end namespace ode
#endif
