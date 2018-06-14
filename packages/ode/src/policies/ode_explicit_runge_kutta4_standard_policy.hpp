
#ifndef ODE_EXPLICIT_RUNGEKUTTA4_STANDARD_POLICY_HPP_
#define ODE_EXPLICIT_RUNGEKUTTA4_STANDARD_POLICY_HPP_

#include "ode_ConfigDefs.hpp"
#include "./base/ode_explicit_policy_base.hpp"

namespace ode{
namespace policy{

template<typename state_type, typename residual_type,
	 typename model_type, typename time_type>
class explicitRungeKutta4StandardResidual
  : public explicitResidualPolicyBase< explicitRungeKutta4StandardResidual<state_type,residual_type,
										      model_type,time_type>,
						  state_type, residual_type,model_type, time_type>
{
 private:
  void computeImpl(const state_type & y, residual_type & R,
		   model_type & model, time_type t){
    model.residual(y, R, t);
  }
};
  
}//end namespace polices
}//end namespace ode  
#endif 
