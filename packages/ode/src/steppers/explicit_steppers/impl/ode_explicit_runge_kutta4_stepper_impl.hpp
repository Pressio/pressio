
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
				 residual_policy_type> >
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

protected:
  using stepper_base_t::model_;
  using stepper_base_t::residual_obj_;
  
protected:
  template < typename T = model_type,
  	     typename U = residual_policy_type,
	     typename... Args>
  explicitRungeKutta4StepperImpl(T & model,
				 U & res_policy_obj,
				 Args&&... rest)
    : stepper_base_t(model, res_policy_obj),
      rhs1_(std::forward<Args>(rest)...),
      rhs2_(std::forward<Args>(rest)...),
      rhs3_(std::forward<Args>(rest)...),
      rhs4_(std::forward<Args>(rest)...),
      y_tmp_(std::forward<Args>(rest)...){}

  explicitRungeKutta4StepperImpl(){};
  ~explicitRungeKutta4StepperImpl() = default;

protected:
  template<typename step_t>
  void doStepImpl(state_type & y,
		  time_type t,
		  time_type dt,
		  step_t step)
  {
    //static const scalar_type val1 = static_cast< scalar_type >( 1 );
    const time_type dt_half = dt / static_cast< scalar_type >(2);
    const time_type t_phalf = t + dt_half;
    const time_type dt6 = dt / static_cast< scalar_type >( 6 );
    const time_type dt3 = dt / static_cast< scalar_type >( 3 );

    auto ySz = sizer_type::getSize(y);
    if(sizer_type::getSize(y_tmp_) == 0)
      sizer_type::resize(y_tmp_, ySz);
    if(sizer_type::getSize(rhs1_) == 0)
      sizer_type::resize(rhs1_, ySz);
    if(sizer_type::getSize(rhs2_) == 0)
      sizer_type::resize(rhs2_, ySz);
    if(sizer_type::getSize(rhs3_) == 0)
      sizer_type::resize(rhs3_, ySz);
    if(sizer_type::getSize(rhs4_) == 0)
      sizer_type::resize(rhs4_, ySz);
    
    // ----------
    // stage 1: 
    // ----------
    // rhs1_(y_n,t)
    residual_obj_->compute(y, rhs1_, *model_, t);
    // y_tmp_ = y_n + rhs1_*dt/2
    for (decltype(ySz) i=0; i<ySz; i++){
      y_tmp_[i] = y[i] + dt_half*rhs1_[i];
    }

    // ----------
    // stage 2: 
    // ----------
    // rhs2_
    residual_obj_->compute(y_tmp_, rhs2_, *model_, t_phalf);
    // y_tmp_ = y_n + rhs2_*dt/2
    for (decltype(ySz) i=0; i<ySz; i++){
      y_tmp_[i] = y[i] + dt_half*rhs2_[i];
    }

    // ----------
    // stage 3: 
    // ----------
    // rhs3_
    residual_obj_->compute(y_tmp_, rhs3_, *model_, t_phalf);
    //y_tmp_ = y_n + rhs3_*dt/2
    for (decltype(ySz) i=0; i<ySz; i++){
      y_tmp_[i] = y[i] + dt*rhs3_[i];
    }

    // ----------
    // stage 4: 
    // ----------
    // rhs4_
    residual_obj_->compute(y_tmp_, rhs4_, *model_, t + dt);
    //x += dt/6 * ( k1 + 2 * k2 + 2 * k3 + k4 )
    for (decltype(ySz) i=0; i < ySz; i++)
    {
      y[i] += dt6*rhs1_[i] + dt3*rhs2_[i] + dt3*rhs3_[i] + dt6*rhs4_[i];
    }

  }
  //----------------------------------------------------------------
  
private:
  friend stepper_base_t;

  residual_type rhs1_;
  residual_type rhs2_;
  residual_type rhs3_;
  residual_type rhs4_;
  state_type y_tmp_;
  // inherited: model_, residual_obj_
  
}; //end class

}//end namespace impl
}//end namespace ode  
#endif
