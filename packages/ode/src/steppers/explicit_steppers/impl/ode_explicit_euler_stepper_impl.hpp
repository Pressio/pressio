
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

  //TODO: sfinae out for special cases, like native operators
  void doStepImpl(state_type & y,
		  ode::details::time_type t,
		  ode::details::time_type dt )
  {
    size_t ySz = sizer_type::getSize(y);

    // rhs has to be same size of lhs
    sizer_type::resize(RHS_, ySz);

    //eval RHS
    this->residual_policy_obj_->compute(y, RHS_,
					*(this->model_),
					t, ySz, ySz);
    
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
