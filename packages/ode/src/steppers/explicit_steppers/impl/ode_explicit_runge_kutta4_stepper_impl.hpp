
#ifndef ODE_EXPLICIT_RUNGEKUTTA4_STEPPER_IMPL_HPP_
#define ODE_EXPLICIT_RUNGEKUTTA4_STEPPER_IMPL_HPP_

#include "../base/ode_explicit_stepper_base.hpp"
#include "../../ode_butcher_tableau.hpp"
#include "ode_runge_kutta_storage.hpp"

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
    private rungeKuttaStorage<state_type, residual_type, 4>
{

private:
  using stepper_t = explicitRungeKutta4StepperImpl<state_type,
						   residual_type,
						   scalar_type,
						   model_type,
						   time_type,
						   sizer_type,
						   residual_policy_type>;
  
  using stepper_base_t = explicitStepperBase<stepper_t>;
  using storage_base_t = rungeKuttaStorage<state_type, residual_type, 4>;
  
protected:
  using stepper_base_t::model_;
  using stepper_base_t::residual_obj_;
  using storage_base_t::rhs_;
  
protected:
  template < typename T = model_type,
  	     typename U = residual_policy_type,
	     typename... Args>
  explicitRungeKutta4StepperImpl(T & model,
				 U & res_policy_obj,
				 Args&&... rest)
    : stepper_base_t(model, res_policy_obj),
      rungeKuttaStorage<state_type,residual_type,4>(std::forward<Args>(rest)...),
    y_tmp_(std::forward<Args>(rest)...){}

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
    if(sizer_type::getSize(y_tmp_) == 0)
      sizer_type::matchSize(y_tmp_, y);

    const time_type dt_half = dt / static_cast< scalar_type >(2);
    const time_type t_phalf = t + dt_half;
    const time_type dt6 = dt / static_cast< scalar_type >( 6 );
    const time_type dt3 = dt / static_cast< scalar_type >( 3 );
    
    if(sizer_type::getSize(y_tmp_) == 0)
      sizer_type::matchSize(y_tmp_, y);
    if(sizer_type::getSize(rhs_[0]) == 0)
      sizer_type::matchSize(rhs_[0], y);
    if(sizer_type::getSize(rhs_[1]) == 0)
      sizer_type::matchSize(rhs_[1], y);
    if(sizer_type::getSize(rhs_[2]) == 0)
      sizer_type::matchSize(rhs_[2], y);
    if(sizer_type::getSize(rhs_[3]) == 0)
      sizer_type::matchSize(rhs_[3], y);
    
    // ----------
    // stage 1: 
    // ----------
    // rhs_[0](y_n,t)
    residual_obj_->compute(y, rhs_[0], *model_, t);
    // y_tmp_ = y_n + rhs_[0]*dt/2
    for (decltype(ySz) i=0; i<ySz; i++){
      y_tmp_[i] = y[i] + dt_half*rhs_[0][i];
    }

    // ----------
    // stage 2: 
    // ----------
    // rhs_[1]
    residual_obj_->compute(y_tmp_, rhs_[1], *model_, t_phalf);
    // y_tmp_ = y_n + rhs_[1]*dt/2
    for (decltype(ySz) i=0; i<ySz; i++){
      y_tmp_[i] = y[i] + dt_half*rhs_[1][i];
    }

    // ----------
    // stage 3: 
    // ----------
    // rhs_[2]
    residual_obj_->compute(y_tmp_, rhs_[2], *model_, t_phalf);
    //y_tmp_ = y_n + rhs_[2]*dt/2
    for (decltype(ySz) i=0; i<ySz; i++){
      y_tmp_[i] = y[i] + dt*rhs_[2][i];
    }

    // ----------
    // stage 4: 
    // ----------
    // rhs_[3]
    residual_obj_->compute(y_tmp_, rhs_[3], *model_, t + dt);
    //x += dt/6 * ( k1 + 2 * k2 + 2 * k3 + k4 )
    for (decltype(ySz) i=0; i < ySz; i++)
    {
      y[i] += dt6*rhs_[0][i] + dt3*rhs_[1][i] + dt3*rhs_[2][i] + dt6*rhs_[3][i];
    }
  }//end doStep

private:
  friend stepper_base_t;

  state_type y_tmp_;
  // inherited: model_, residual_obj_
  
}; //end class

}//end namespace impl
}//end namespace ode  
#endif
