
#ifndef ODE_STEPPERS_EXPLICIT_STEPPERS_IMPL_EXPLICIT_EULER_STEPPER_IMPL_HPP_
#define ODE_STEPPERS_EXPLICIT_STEPPERS_IMPL_EXPLICIT_EULER_STEPPER_IMPL_HPP_

#include "../base/ode_explicit_stepper_base.hpp"

namespace rompp{
namespace ode{
namespace impl{
    
template<typename ode_state_type,
	 typename model_type,	
	 typename ode_residual_type,
	 typename residual_policy_type
	 >
class ExplicitEulerStepperImpl<ode_state_type,
			       model_type,
			       ode_residual_type,
			       residual_policy_type>
  : public ExplicitStepperBase<
  ExplicitEulerStepperImpl<ode_state_type,	
			   model_type,
			   ode_residual_type,
			   residual_policy_type>>,
    private OdeStorage<ode_state_type, ode_residual_type, 0, 1>,
    private ExpOdeAuxData<model_type, residual_policy_type>{

  static_assert( meta::is_legitimate_explicit_residual_policy<
		 residual_policy_type>::value ||
		 meta::is_explicit_euler_residual_standard_policy<
		 residual_policy_type>::value,
"EXPLICIT EULER RESIDUAL_POLICY NOT ADMISSIBLE, \
MAYBE NOT A CHILD OF ITS BASE OR DERIVING FROM WRONG BASE");

  using this_t = ExplicitEulerStepperImpl<
    ode_state_type, model_type, ode_residual_type, residual_policy_type>;
  using stepper_base_t = ExplicitStepperBase<this_t>;
  using storage_base_t = OdeStorage<ode_state_type, ode_residual_type, 0, 1>;
  using auxdata_base_t = ExpOdeAuxData<model_type, residual_policy_type>;
  using scalar_type  = typename core::details::traits<ode_state_type>::scalar_t;
  using scalar_t2  = typename core::details::traits<ode_residual_type>::scalar_t;
  static_assert(std::is_same<scalar_type, scalar_t2>::value,
		"Not maching scalar types");
    
protected:
  using storage_base_t::auxRHS_;
  using auxdata_base_t::model_;
  using auxdata_base_t::residual_obj_;

protected:

  template <typename T1 = model_type,
  	    typename T2 = residual_policy_type,
	    typename T3 = ode_state_type,
	    typename T4 = ode_residual_type,
	    typename... Args>
  ExplicitEulerStepperImpl(const T1 & model,
			   const T2 & res_policy_obj,
			   const T3 & y0,
			   const T4 & r0,
			   Args&&... rest)
    : storage_base_t(r0 /*,std::forward<Args>(rest)...*/),
      auxdata_base_t(model, res_policy_obj){}

  template <typename T1 = residual_policy_type,
	    typename T2 = ode_state_type,
	    typename T3 = ode_residual_type,
	    typename... Args>
  ExplicitEulerStepperImpl(const T1 & res_policy_obj,
			   const T2 & y0,
			   const T3 & r0,
			   Args&&... rest)
    : storage_base_t(r0 /*,std::forward<Args>(rest)...*/),
      auxdata_base_t(res_policy_obj){}
  
  ExplicitEulerStepperImpl() = delete;
  ~ExplicitEulerStepperImpl() = default;

public:
  
  template<typename step_t>
  void operator()(ode_state_type & y, scalar_type t,
		  scalar_type dt, step_t step){
    
    if ( auxRHS_[0].empty() )
      auxRHS_[0].matchLayoutWith(y);

    //eval RHS
    this->evalRHS(y, auxRHS_[0], t);
    
    // y = y + dt * rhs
    y += dt * auxRHS_[0];
  }
  //-------------------------------------------------------

protected:
  
  template<typename T = model_type,
	   typename std::enable_if<
	     std::is_void<T>::value
	     >::type * = nullptr>
  void evalRHS(const ode_state_type & y,
	       ode_residual_type & RHS,
	       scalar_type t){
    (*residual_obj_)(y, RHS, t);
  }
  //-------------------------------------------------------
  
  template<typename T = model_type,
	   typename std::enable_if<
	     !std::is_void<T>::value
	     >::type * = nullptr>
  void evalRHS(const ode_state_type & y,
	       ode_residual_type & RHS,
	       scalar_type t){
    (*residual_obj_)(y, RHS, *model_, t);
  }
  //----------------------------------------------------------------
  
private:
  friend stepper_base_t;
  
}; //end class

}//end namespace impl
}//end namespace ode  
}//end namespace rompp
#endif
