
#ifndef ODE_EXPLICIT_RUNGEKUTTA4_STEPPER_IMPL_HPP_
#define ODE_EXPLICIT_RUNGEKUTTA4_STEPPER_IMPL_HPP_

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
class explicitRungeKutta4StepperImpl<state_type,
				     residual_type,
				     scalar_type,
				     model_type,
				     time_type,
				     sizer_type,
				     residual_policy_type>
  : public explicitStepperBase<
  explicitRungeKutta4StepperImpl<state_type,
				 residual_type,
				 scalar_type,
				 model_type,
				 time_type,
				 sizer_type,
				 residual_policy_type> >,
    private odeStorage<state_type, residual_type, 1, 4>
{
  static_assert( meta::isLegitimateExplicitResidualPolicy<
		 residual_policy_type>::value ||
		 meta::isExplicitRungeKutta4ResidualStandardPolicy<
		 residual_policy_type>::value,
     "EXPLICIT RUNGEKUTTA4 RESIDUAL_POLICY NOT ADMISSIBLE, \
MAYBE NOT A CHILD OF ITS BASE OR DERIVING FROM WRONG BASE");

private:
  using stepper_t = explicitRungeKutta4StepperImpl<
  state_type, residual_type, scalar_type,
  model_type, time_type, sizer_type, residual_policy_type>;
  
  using stepper_base_t = explicitStepperBase<stepper_t>;
  using storage_base_t = odeStorage<state_type, residual_type, 1, 4>;
  
protected:
  using stepper_base_t::model_;
  using stepper_base_t::residual_obj_;
  using storage_base_t::auxStates_;
  using storage_base_t::auxRHS_;
  
protected:
  template < typename T = model_type,
  	     typename U = residual_policy_type,
	     typename... Args>
  explicitRungeKutta4StepperImpl(T & model,
				 U & res_policy_obj,
				 Args&&... rest)
    : stepper_base_t(model, res_policy_obj),
      storage_base_t(std::forward<Args>(rest)...){}

  explicitRungeKutta4StepperImpl() = delete;
  ~explicitRungeKutta4StepperImpl() = default;

protected:
  template<typename step_t>
  void doStepImpl(state_type & y,
		  time_type t,
		  time_type dt,
		  step_t step)
  {
    auto ySz = sizer_type::getSize(y);
    if(sizer_type::getSize(auxStates_[0]) == 0)
      sizer_type::matchSize(y, auxStates_[0]);

    const time_type dt_half = dt / static_cast< scalar_type >(2);
    const time_type t_phalf = t + dt_half;
    const time_type dt6 = dt / static_cast< scalar_type >( 6 );
    const time_type dt3 = dt / static_cast< scalar_type >( 3 );

    auto & ytmp = auxStates_[0];
    
    if(sizer_type::getSize(ytmp) == 0)
      sizer_type::matchSize(y, ytmp);
    if(sizer_type::getSize(auxRHS_[0]) == 0)
      sizer_type::matchSize(y, auxRHS_[0]);
    if(sizer_type::getSize(auxRHS_[1]) == 0)
      sizer_type::matchSize(y, auxRHS_[1]);
    if(sizer_type::getSize(auxRHS_[2]) == 0)
      sizer_type::matchSize(y, auxRHS_[2]);
    if(sizer_type::getSize(auxRHS_[3]) == 0)
      sizer_type::matchSize(y, auxRHS_[3]);
    
    // ----------
    // stage 1: 
    // ----------
    // rhs_[0](y_n,t)
    residual_obj_->compute(y, auxRHS_[0], *model_, t);
    // ytmp = y_n + auxRHS_[0]*dt/2
    for (decltype(ySz) i=0; i<ySz; i++){
      ytmp[i] = y[i] + dt_half*auxRHS_[0][i];
    }

    // ----------
    // stage 2: 
    // ----------
    // rhs_[1]
    residual_obj_->compute(ytmp, auxRHS_[1], *model_, t_phalf);
    // ytmp = y_n + auxRHS_[1]*dt/2
    for (decltype(ySz) i=0; i<ySz; i++){
      ytmp[i] = y[i] + dt_half*auxRHS_[1][i];
    }

    // ----------
    // stage 3: 
    // ----------
    // auxRHS_[2]
    residual_obj_->compute(ytmp, auxRHS_[2], *model_, t_phalf);
    //ytmp = y_n + auxRHS_[2]*dt/2
    for (decltype(ySz) i=0; i<ySz; i++){
      ytmp[i] = y[i] + dt*auxRHS_[2][i];
    }

    // ----------
    // stage 4: 
    // ----------
    // auxRHS_[3]
    residual_obj_->compute(ytmp, auxRHS_[3], *model_, t + dt);
    //x += dt/6 * ( k1 + 2 * k2 + 2 * k3 + k4 )
    for (decltype(ySz) i=0; i < ySz; i++)
    {
      y[i] += dt6*auxRHS_[0][i] + dt3*auxRHS_[1][i] +
	dt3*auxRHS_[2][i] + dt6*auxRHS_[3][i];
    }
  }//end doStep

private:
  friend stepper_base_t;
  // inherited: model_, residual_obj_
  
}; //end class

}//end namespace impl
}//end namespace ode  
#endif
