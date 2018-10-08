
#ifndef ODE_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_HPP_
#define ODE_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_HPP_

#include "ode_implicit_stepper_helper_info.hpp"

namespace rompp{ namespace ode{ 
      
template<ImplicitSteppersEnum whichone, typename... Args>
class ImplicitStepper
  : public impl::implicit_stepper_helper_info<whichone,
					      Args...>::base_impl_type{

  using info_t = impl::implicit_stepper_helper_info<whichone, Args...>;
  using state_t = typename info_t::state_type;
  using model_t = typename info_t::model_type;
  using res_t = typename info_t::res_type;
  using jac_t = typename info_t::jac_type;
  using residual_pol_t = typename info_t::residual_policy_type;
  using residual_pol_std_t = typename info_t::res_std_pol_type;
  using jacobian_pol_t = typename info_t::jacobian_policy_type;
  using jacobian_pol_std_t = typename info_t::jac_std_pol_type;
  using auxiliary_stepper_t = typename info_t::auxiliary_stepper_type;
  using base_impl_t = typename info_t::base_impl_type;

//this needs to be public, it is detected by integrators
public:
  using base_t = base_impl_t;
  
public:

  // some constructors need y0, and some need y0/r0 
  // because some multi-stage methods need history of RHS
  // e.g. backward Euler only needs y, y_n-1
  // e.g. BDF2 only needs y, y_n-1, y_n-2

  // passing: model, initial state
  // policy is standard for residual and jacobian
  ImplicitStepper(const model_t & model,
		  const state_t & y0)
    : base_impl_t(model,
		  res_policy_obj_,
		  jac_policy_obj_,
		  y0){}

  // passing: model, initial state, aux_stepper
  // policy is standard for residual and jacobian
  template <typename T = auxiliary_stepper_t,
  	    core::meta::enable_if_t<
  	      not std::is_void<T>::value
  	      > * = nullptr>
  ImplicitStepper(const model_t & model,
  		  const state_t & y0,
  		  T & auxStObj)
    : base_impl_t(model,
  		  res_policy_obj_,
  		  jac_policy_obj_,
  		  auxStObj, y0){}
  
  // passing: model, initial state, and policies. 
  ImplicitStepper(const model_t & model,
		  const residual_pol_t & resPolicyObj,
		  const jacobian_pol_t & jacPolicyObj,
		  const state_t & y0)
    : base_impl_t(model,
      resPolicyObj,
      jacPolicyObj,
      y0){}

  // passing: model, initial state, initial residual
  // policy is standard for residual and jacobian
  ImplicitStepper(const model_t & model,
		  const state_t & y0,
		  const res_t & r0)
    : base_impl_t(model,
      res_policy_obj_,
      jac_policy_obj_,
      y0, r0){}

  // passing: model, initial state, 
  // initial residual and policies. 
  ImplicitStepper(const model_t & model,
  		  const residual_pol_t & resPolicyObj,
		  const jacobian_pol_t & jacPolicyObj,
  		  const state_t & y0,
  		  const res_t & r0)
    : base_impl_t(model,
		  resPolicyObj,
		  jacPolicyObj,
		  y0, r0){}
  
  ImplicitStepper() = delete;
  virtual ~ImplicitStepper() = default;

private:
  // not used if policy is passed from outside
  residual_pol_t res_policy_obj_;
  jacobian_pol_t jac_policy_obj_;

};//end class

}} // end namespace rompp::ode
#endif 
