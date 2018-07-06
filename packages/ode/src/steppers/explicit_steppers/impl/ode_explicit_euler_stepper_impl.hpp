
#ifndef ODE_EXPLICIT_EULER_STEPPER_IMPL_HPP_
#define ODE_EXPLICIT_EULER_STEPPER_IMPL_HPP_

#include "../base/ode_explicit_stepper_base.hpp"

namespace ode{
namespace impl{
    
template<typename state_type,
	 typename residual_type,
	 typename scalar_type,
	 typename model_type,	
	 typename time_type,
	 typename sizer_type,
	 typename residual_policy_type
	 >
class explicitEulerStepperImpl<state_type,
			       residual_type,
			       scalar_type,
			       model_type,
			       time_type,
			       sizer_type,
			       residual_policy_type>
  : public explicitStepperBase<
  explicitEulerStepperImpl<state_type,
			   residual_type,
			   scalar_type,
			   model_type,
			   time_type,
			   sizer_type,
			   residual_policy_type> >
{  
private:
  using stepper_t = explicitEulerStepperImpl<state_type,
					     residual_type,
					     scalar_type,
					     model_type,
					     time_type,
					     sizer_type,
					     residual_policy_type>;
  using stepper_base_t = explicitStepperBase<stepper_t>;

protected:
  using stepper_base_t::model_;
  using stepper_base_t::residual_obj_;
  using stepper_base_t::mysizer_;
  
protected:
  template < typename T = model_type,
  	     typename U = residual_policy_type,
	     typename... Args>
  explicitEulerStepperImpl(T & model,
			   U & res_policy_obj,
			   Args&&... rest)
    : stepper_base_t(model, res_policy_obj),
      RHS_(std::forward<Args>(rest)...){}

  explicitEulerStepperImpl() = delete;
  ~explicitEulerStepperImpl() = default;

protected:
  template<typename step_t>
  void doStepImpl(state_type & y,
		  time_type t,
		  time_type dt,
		  step_t step)
  {
    size_t ySz = sizer_type::getSize(y);

    //eval RHS
    residual_obj_->compute(y, RHS_, *model_, t, mysizer_);
    
    //out = in + dt * rhs
    for (size_t i=0; i < ySz; i++){
      y[i] += dt*RHS_[i];
    }
  }
  //----------------------------------------------------------------
  
private:
  friend stepper_base_t;
  residual_type RHS_;
  // inherited: model_, residual_policy_obj_
  
}; //end class

}//end namespace impl
}//end namespace ode  
#endif
