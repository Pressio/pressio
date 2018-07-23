
#ifndef ODE_STEPPERS_IMPLICIT_STEPPERS_IMPL_IMPLICIT_EULER_STEPPER_IMPL_HPP_
#define ODE_STEPPERS_IMPLICIT_STEPPERS_IMPL_IMPLICIT_EULER_STEPPER_IMPL_HPP_

#include "../base/ode_implicit_stepper_base.hpp"
#include "../../../policies/meta/ode_implicit_euler_policies_meta.hpp"

namespace ode{
namespace impl{

template<typename state_type,
	 typename residual_type,
	 typename jacobian_type,
	 typename scalar_type,
	 typename model_type,
	 typename sizer_type,
	 typename residual_policy_type,
	 typename jacobian_policy_type>
class ImplicitEulerStepperImpl<state_type,
			       residual_type,
			       jacobian_type,
			       scalar_type,
			       model_type,
			       sizer_type,
			       residual_policy_type,
			       jacobian_policy_type>
  : public ImplicitStepperBase<
	    ImplicitEulerStepperImpl<state_type, residual_type,
				     jacobian_type, scalar_type,
				     model_type, sizer_type,
				     residual_policy_type,
				     jacobian_policy_type> >,
    private OdeStorage<state_type, residual_type, 1>,
    private ImpOdeAuxData<model_type, scalar_type,
			  residual_policy_type,
			  jacobian_policy_type>  
{

  static_assert( meta::is_legitimate_implicit_euler_residual_policy<
		 residual_policy_type>::value ||
		 meta::is_implicit_euler_residual_standard_policy<
		 residual_policy_type>::value,
		 "IMPLICIT EULER RESIDUAL_POLICY NOT ADMISSIBLE, \
MAYBE NOT A CHILD OF ITS BASE OR DERIVING FROM WRONG BASE");

  static_assert( meta::is_legitimate_implicit_euler_jacobian_policy<
		 jacobian_policy_type>::value,
		 "IMPLICIT EULER JACOBIAN_POLICY NOT ADMISSIBLE, \
MAYBE NOT A CHILD OF ITS BASE OR DERIVING FROM WRONG BASE");
  
  using stepper_t = ImplicitEulerStepperImpl<state_type,
					     residual_type,
					     jacobian_type,
					     scalar_type,
					     model_type,
					     sizer_type,
					     residual_policy_type,
					     jacobian_policy_type>;
  using stepper_base_t = ImplicitStepperBase<stepper_t>;
  using storage_base_t = OdeStorage<state_type, residual_type, 1>;
  using auxdata_base_t = ImpOdeAuxData<model_type, scalar_type,
				       residual_policy_type,
				       jacobian_policy_type>;

protected:
  using storage_base_t::auxStates_;

  using auxdata_base_t::model_;
  using auxdata_base_t::residual_obj_;
  using auxdata_base_t::jacobian_obj_;
  using auxdata_base_t::t_;
  using auxdata_base_t::dt_;
  
protected:
  template < typename M = model_type,
	     typename U = residual_policy_type,
	     typename T = jacobian_policy_type,
	     typename... Args>
  ImplicitEulerStepperImpl(M & model,
			   U & res_policy_obj,
			   T & jac_policy_obj,
			   Args&& ... rest)
    : storage_base_t(std::forward<Args>(rest)...),
      auxdata_base_t(model, res_policy_obj, jac_policy_obj)
  {}
    
  ImplicitEulerStepperImpl() = delete;
  ~ImplicitEulerStepperImpl() = default;
  
protected:
  template<typename solver_type, typename step_t>
  void doStepImpl(state_type & y, scalar_type t,
		  scalar_type dt, step_t step, solver_type & solver)
  {
    dt_ = dt;
    t_ = t;

    // store previous state = y;
    auxStates_[0] = y;
    solver.solve(y, *this);

  }//end doStepImpl

  void residualImpl(const state_type & y, residual_type & R){
    residual_obj_->compute(y, R, auxStates_, *model_, t_, dt_);
  }

  void jacobianImpl(const state_type & y, jacobian_type & J)
  {
    jacobian_obj_->compute(y, J, *model_, t_, dt_);
  }
  
private:
  friend stepper_base_t;

}; //end class

}//end namespace impl
}//end namespace ode  
#endif 
