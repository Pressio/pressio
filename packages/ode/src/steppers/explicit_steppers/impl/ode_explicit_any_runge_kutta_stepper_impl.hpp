
#ifndef ODE_EXPLICIT_ANY_RUNGEKUTTA4_STEPPER_IMPL_HPP_
#define ODE_EXPLICIT_ANY_RUNGEKUTTA4_STEPPER_IMPL_HPP_

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
	 typename residual_policy_type,
	 typename butcher_table_type>
class explicitAnyRungeKuttaStepperImpl<state_type,
				       residual_type,
				       scalar_type,
				       model_type,
				       time_type,
				       sizer_type,
				       residual_policy_type,
				       butcher_table_type>
  : public explicitStepperBase<
  explicitAnyRungeKuttaStepperImpl<state_type,
				   residual_type,
				   scalar_type,
				   model_type,
				   time_type,
				   sizer_type,
				   residual_policy_type,
				   butcher_table_type> >,
    private rungeKuttaStorage<state_type, residual_type, 4>
{

private:
  using stepper_t = explicitAnyRungeKuttaStepperImpl<state_type,
						   residual_type,
						   scalar_type,
						   model_type,
						   time_type,
						   sizer_type,
						   residual_policy_type,
						   butcher_table_type>;
  
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
  explicitAnyRungeKuttaStepperImpl(T & model,
				 U & res_policy_obj,
				 butcher_table_type bT,
				 Args&&... rest)
    : stepper_base_t(model, res_policy_obj),
      storage_base_t(std::forward<Args>(rest)...),
      y_tmp_(std::forward<Args>(rest)...),
      bT_(bT){}

  explicitAnyRungeKuttaStepperImpl() = delete;
  ~explicitAnyRungeKuttaStepperImpl() = default;

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

    int nStages = 4;
    for (int iStage=1; iStage<=nStages; iStage++)
    {
      if(sizer_type::getSize(rhs_[iStage-1]) == 0)
    	sizer_type::matchSize(rhs_[iStage-1], y);

      time_type tnow = t + dt*bT_.c(iStage);	
      if (iStage==1)
    	residual_obj_->compute(y, rhs_[iStage-1], *model_, tnow);
      else
    	residual_obj_->compute(y_tmp_, rhs_[iStage-1], *model_, tnow);

      if (iStage<nStages){
    	for (decltype(ySz) i=0; i<ySz; i++){
    	  auto sum = static_cast<scalar_type>(0);
    	  for (int j=1; j<=iStage; j++){
    	    sum += dt*bT_.a(iStage+1,j)*rhs_[j-1][i];
    	  }
    	  y_tmp_[i] = y[i] + sum;
    	}
      }
    }//end stages for

    for (decltype(ySz) i=0; i<ySz; i++){
      for (int j=1; j<=nStages; j++)
    	y[i] += dt*bT_.b(j)*rhs_[j-1][i];
    }
  }
  
private:
  friend stepper_base_t;

  state_type y_tmp_;
  butcher_table_type bT_;
  // inherited: model_, residual_obj_
  
}; //end class

}//end namespace impl
}//end namespace ode  
#endif
