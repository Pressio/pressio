
#ifndef ODE_EXPLICIT_EULER_STEPPER_IMPL_HPP_
#define ODE_EXPLICIT_EULER_STEPPER_IMPL_HPP_

#include "../base/ode_explicit_stepper_base.hpp"
#include "../../../policies/meta/ode_explicit_euler_policies_meta.hpp"

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
class ExplicitEulerStepperImpl<state_type,
			       residual_type,
			       scalar_type,
			       model_type,
			       time_type,
			       sizer_type,
			       residual_policy_type>
  : public ExplicitStepperBase<
  ExplicitEulerStepperImpl<state_type,
			   residual_type,
			   scalar_type,
			   model_type,
			   time_type,
			   sizer_type,
			   residual_policy_type> >,
    private OdeStorage<state_type, residual_type, 0, 1>
{  
  static_assert( meta::is_legitimate_explicit_residual_policy<
		 residual_policy_type>::value ||
		 meta::is_explicit_euler_residual_standard_policy<
		 residual_policy_type>::value,
	  "EXPLICIT EULER RESIDUAL_POLICY NOT ADMISSIBLE, \
MAYBE NOT A CHILD OF ITS BASE OR DERIVING FROM WRONG BASE");

  using stepper_t = ExplicitEulerStepperImpl<
  state_type, residual_type, scalar_type,
  model_type, time_type, sizer_type, residual_policy_type>;

  using stepper_base_t = ExplicitStepperBase<stepper_t>;
  using storage_base_t = OdeStorage<state_type, residual_type, 0, 1>;

protected:
  using stepper_base_t::model_;
  using stepper_base_t::residual_obj_;
  using storage_base_t::auxRHS_;
  
protected:
  template < typename T = model_type,
  	     typename U = residual_policy_type,
	     typename... Args>
  ExplicitEulerStepperImpl(T & model,
			   U & res_policy_obj,
			   Args&&... rest)
    : stepper_base_t(model, res_policy_obj),
      storage_base_t(std::forward<Args>(rest)...){}

  ExplicitEulerStepperImpl() = delete;
  ~ExplicitEulerStepperImpl() = default;

protected:
  template<typename step_t>
  void doStepImpl(state_type & y, time_type t,
		  time_type dt, step_t step)
  {
    auto ySz = sizer_type::getSize(y);
    if (sizer_type::getSize(auxRHS_[0]) == 0)
      sizer_type::matchSize(y, auxRHS_[0]);

    //eval RHS
    residual_obj_->compute(y, auxRHS_[0], *model_, t);
    
    // //out = in + dt * rhs
    for (decltype(ySz) i=0; i < ySz; i++){
      y[i] += dt*auxRHS_[0][i];
    }
  }
  //----------------------------------------------------------------
  
private:
  friend stepper_base_t;
  // inherited: model_, residual_obj_
  
}; //end class

}//end namespace impl
}//end namespace ode  
#endif
