
#ifndef ODE_IMPLICIT_EULER_STEPPER_IMPL_HPP_
#define ODE_IMPLICIT_EULER_STEPPER_IMPL_HPP_

#include "../base/ode_implicit_stepper_base.hpp"

namespace ode{
namespace impl{

template<typename state_type,
	 typename residual_type,
	 typename jacobian_type,
	 typename scalar_type,
	 typename model_type,
	 typename time_type,
	 typename sizer_type,
	 typename solver_policy_type,
	 typename residual_policy_type,
	 typename jacobian_policy_type>
class implicitEulerStepperImpl<state_type,
			       residual_type,
			       jacobian_type,
			       scalar_type,
			       model_type,
			       time_type,
			       sizer_type,
			       solver_policy_type,
			       residual_policy_type,
			       jacobian_policy_type>
  : public implicitStepperBase<
  implicitEulerStepperImpl<state_type, residual_type,
			   jacobian_type, scalar_type,
			   model_type, time_type,
			   sizer_type, solver_policy_type,
			   residual_policy_type,
			   jacobian_policy_type> >
{
  static_assert( meta::isLegitimateImplicitEulerResidualPolicy<residual_policy_type>::value,
		 "IMPLICIT EULER RESIDUAL_POLICY NOT ADMISSIBLE, MAYBE NOT A CHILD OF ITS BASE OR DERIVING FROM WRONG BASE");
  static_assert( meta::isLegitimateImplicitEulerJacobianPolicy<jacobian_policy_type>::value,
		 "IMPLICIT EULER JACOBIAN_POLICY NOT ADMISSIBLE, MAYBE NOT A CHILD OF ITS BASE OR DERIVING FROM WRONG BASE");
  
private:
  using stepper_t = implicitEulerStepperImpl<state_type,
					     residual_type,
					     jacobian_type,
					     scalar_type,
					     model_type,
					     time_type,
					     sizer_type,
					     solver_policy_type,
					     residual_policy_type,
					     jacobian_policy_type>;
  using stepper_base_t = implicitStepperBase<stepper_t>;

protected:
  template < typename M = model_type,
	     typename S = solver_policy_type,
	     typename U = residual_policy_type,
	     typename T = jacobian_policy_type,
	     typename... Args>
  implicitEulerStepperImpl(M & model,
			   S & solver,
			   U & res_policy_obj,
			   T & jac_policy_obj,
			   Args&& ... rest)
    : stepper_base_t(model, res_policy_obj, jac_policy_obj), solver_(&solver),
      y_nm1_(std::forward<Args>(rest)...)
  {}
    
  implicitEulerStepperImpl() = delete;
  ~implicitEulerStepperImpl() = default;
  
protected:
  template<typename step_t>
  void doStepImpl(state_type & y, time_type t, time_type dt, step_t step )
  {
    y_nm1_ = y;
    dt_ = dt;
    t_ = t;
    solver_->solve(y, *this);
  }//end doStepImpl

  void residualImpl(const state_type & y, state_type & R){
    this->residual_policy_obj_->compute(y, y_nm1_, R,
					*(this->model_), t_, dt_);
  }
  void jacobianImpl(const state_type & y, jacobian_type & J){
    this->jacobian_policy_obj_->compute(y, J, *(this->model_), t_, dt_);
  }

private:
  friend stepper_base_t;
  solver_policy_type * solver_;
  time_type t_;
  time_type dt_;
  state_type y_nm1_;

}; //end class

}//end namespace impl
}//end namespace ode  
#endif 
