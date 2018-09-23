
#ifndef ODE_STEPPERS_EXPLICIT_STEPPERS_IMPL_EXPLICIT_RUNGEKUTTA4_STEPPER_IMPL_HPP_
#define ODE_STEPPERS_EXPLICIT_STEPPERS_IMPL_EXPLICIT_RUNGEKUTTA4_STEPPER_IMPL_HPP_

#include "../base/ode_explicit_stepper_base.hpp"
#include "../../../policies/meta/ode_explicit_runge_kutta4_policies_meta.hpp"

namespace rompp{
namespace ode{
namespace impl{
  
template<typename state_type,
	 typename ode_residual_type,
	 typename scalar_type,
	 typename model_type,	
	 typename residual_policy_type
	 >
class ExplicitRungeKutta4StepperImpl<state_type,
				     ode_residual_type,
				     scalar_type,
				     model_type,
				     residual_policy_type>
  : public ExplicitStepperBase<
  ExplicitRungeKutta4StepperImpl<state_type,
				 ode_residual_type,
				 scalar_type,
				 model_type,
				 residual_policy_type> >,
    private OdeStorage<state_type, ode_residual_type, 1, 4>,
    private ExpOdeAuxData<model_type, residual_policy_type>
{

  static_assert( meta::is_legitimate_explicit_residual_policy<
		 residual_policy_type>::value ||
		 meta::is_explicit_runge_kutta4_residual_standard_policy<
		 residual_policy_type>::value,
"EXPLICIT RUNGEKUTTA4 RESIDUAL_POLICY NOT ADMISSIBLE, \
MAYBE NOT A CHILD OF ITS BASE OR DERIVING FROM WRONG BASE");

  using stepper_t = ExplicitRungeKutta4StepperImpl<
	   state_type, ode_residual_type, scalar_type,
	   model_type, residual_policy_type>;
  using stepper_base_t = ExplicitStepperBase<stepper_t>;
  using storage_base_t = OdeStorage<state_type, ode_residual_type, 1, 4>;
  using auxdata_base_t = ExpOdeAuxData<model_type, residual_policy_type>;

  using add_op_t = std::plus<scalar_type>;
  
protected:
  using storage_base_t::auxStates_;
  using storage_base_t::auxRHS_;
  using auxdata_base_t::model_;
  using auxdata_base_t::residual_obj_;
  
protected:
  template < typename T1 = model_type,
  	     typename T2 = residual_policy_type,
	     typename T3 = state_type,
	     typename T4 = ode_residual_type,
	     typename... Args>
  ExplicitRungeKutta4StepperImpl(T1 & model,
				 T2 & res_policy_obj,
				 T3 const & y0,
				 T4 const & r0,
				 Args&&... rest)
    : storage_base_t(y0, r0 /*,std::forward<Args>(rest)...*/),
      auxdata_base_t(model, res_policy_obj){}
  
  ExplicitRungeKutta4StepperImpl() = delete;
  ~ExplicitRungeKutta4StepperImpl() = default;

protected:

  template<typename step_t>
  void doStepImpl(state_type & y, scalar_type t,
		  scalar_type dt, step_t step){
    
    auto & ytmp = auxStates_[0];
    if ( ytmp.empty() )
      ytmp.matchLayoutWith(y);
    if(auxRHS_[0].empty() == 0)
      auxRHS_[0].matchLayoutWith(y);
    if(auxRHS_[1].empty() == 0)
      auxRHS_[1].matchLayoutWith(y);
    if(auxRHS_[2].empty() == 0)
      auxRHS_[2].matchLayoutWith(y);
    if(auxRHS_[3].empty() == 0)
      auxRHS_[3].matchLayoutWith(y);

    const scalar_type dt_half = dt / static_cast< scalar_type >(2);
    const scalar_type t_phalf = t + dt_half;
    const scalar_type dt6 = dt / static_cast< scalar_type >( 6 );
    const scalar_type dt3 = dt / static_cast< scalar_type >( 3 );
    
    // ----------
    // stage 1: 
    // ----------
    residual_obj_->compute(y, auxRHS_[0], *model_, t);
    // ytmp = y_n + auxRHS_[0]*dt/2
    ytmp.template inPlaceOp<add_op_t>(1.0,  y, dt_half, auxRHS_[0]);
    
    // ----------
    // stage 2: 
    // ----------
    residual_obj_->compute(ytmp, auxRHS_[1], *model_, t_phalf);
    ytmp.template inPlaceOp<add_op_t>(1.0, y, dt_half, auxRHS_[1]);
    
    // ----------
    // stage 3: 
    // ----------
    residual_obj_->compute(ytmp, auxRHS_[2], *model_, t_phalf);
    //ytmp = y_n + auxRHS_[2]*dt/2
    ytmp.template inPlaceOp<add_op_t>(1.0, y, dt, auxRHS_[2]);

    // ----------
    // stage 4: 
    // ----------
    residual_obj_->compute(ytmp, auxRHS_[3], *model_, t + dt);
    //x += dt/6 * ( k1 + 2 * k2 + 2 * k3 + k4 )
    y.template inPlaceOp<add_op_t>(1.0, dt6, auxRHS_[0],
				   dt3, auxRHS_[1],
				   dt3, auxRHS_[2],
				   dt6, auxRHS_[3]);
    
  }//end doStep

private:
  friend stepper_base_t;
  
}; //end class

}//end namespace impl
}//end namespace ode  
}//end namespace rompp
#endif
