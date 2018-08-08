
#ifndef ODE_STEPPERS_EXPLICIT_STEPPERS_IMPL_EXPLICIT_EULER_STEPPER_IMPL_HPP_
#define ODE_STEPPERS_EXPLICIT_STEPPERS_IMPL_EXPLICIT_EULER_STEPPER_IMPL_HPP_

#include "../base/ode_explicit_stepper_base.hpp"
#include "../../../policies/meta/ode_explicit_euler_policies_meta.hpp"

namespace ode{
namespace impl{
    
template<typename state_type,
	 typename space_residual_type,
	 typename scalar_type,
	 typename model_type,	
	 typename residual_policy_type
	 >
class ExplicitEulerStepperImpl<state_type,
			       space_residual_type,
			       scalar_type,
			       model_type,
			       residual_policy_type>
  : public ExplicitStepperBase<
  ExplicitEulerStepperImpl<state_type,
			   space_residual_type,
			   scalar_type,
			   model_type,
			   residual_policy_type> >,
    private OdeStorage<state_type, space_residual_type, 0, 1>,
    private ExpOdeAuxData<model_type, residual_policy_type>
{  

  static_assert( meta::is_legitimate_explicit_residual_policy<
		 residual_policy_type>::value ||
		 meta::is_explicit_euler_residual_standard_policy<
		 residual_policy_type>::value,
	  "EXPLICIT EULER RESIDUAL_POLICY NOT ADMISSIBLE, \
MAYBE NOT A CHILD OF ITS BASE OR DERIVING FROM WRONG BASE");

  using stepper_t = ExplicitEulerStepperImpl<
    state_type, space_residual_type, scalar_type,
    model_type, residual_policy_type>;

  using stepper_base_t = ExplicitStepperBase<stepper_t>;
  using storage_base_t = OdeStorage<state_type, space_residual_type, 0, 1>;
  using auxdata_base_t = ExpOdeAuxData<model_type, residual_policy_type>;

protected:
  using storage_base_t::auxRHS_;
  using auxdata_base_t::model_;
  using auxdata_base_t::residual_obj_;
  
protected:
  template <typename T1 = model_type,
  	    typename T2 = residual_policy_type,
	    typename T3 = state_type,
	    typename T4 = space_residual_type,
	    typename... Args>
  ExplicitEulerStepperImpl(T1 & model, T2 & res_policy_obj,
			   T3 const & y0, T4 const & r0,
			   Args&&... rest)
    : storage_base_t(r0 /*,std::forward<Args>(rest)...*/),
      auxdata_base_t(model, res_policy_obj){}

  ExplicitEulerStepperImpl() = delete;
  ~ExplicitEulerStepperImpl() = default;

protected:
  template<typename step_t>
  void doStepImpl(state_type & y, scalar_type t,
		  scalar_type dt, step_t step)
  {
    if ( auxRHS_[0].empty() )
      auxRHS_[0].matchLayoutWith(y);

    //eval RHS
    residual_obj_->compute(y, auxRHS_[0], *model_, t);
    
    // y = y + dt * rhs
    y.template inPlaceOp<std::plus<double> >(static_cast<scalar_type>(1.0),
					     dt, auxRHS_[0]);
  }
  //----------------------------------------------------------------
  
private:
  friend stepper_base_t;
  
}; //end class

}//end namespace impl
}//end namespace ode  
#endif
